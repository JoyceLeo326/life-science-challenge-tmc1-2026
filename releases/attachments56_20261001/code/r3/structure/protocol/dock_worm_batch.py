"""Independent R3 7USX experimental local-complete worm condition; original R2 packages remain read-only.

Example: old_env_python dock_batch.py --input <csv> --batch pilot_design --limit 2
Writes raw logs, ligand SDF/PDBQT and poses under --scratch, plus a compact ledger.
Independent pose QC remains a separate review by 05; local checks are provisional.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from portable_config import required_path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from dimorphite_dl import protonate_smiles
from meeko import MoleculePreparation, PDBQTWriterLegacy
import numpy as np
from scipy.spatial import cKDTree
from chemistry_rules import exchangeable_n_indices, stable_identity_smiles, unspecified_potential_stereo, RULE_VERSION

FIRST = required_path("TMC1_R1_ROOT")
DEFAULT_SCRATCH = required_path("TMC1_WORK_ROOT") / "03_compute"
VINA = None  # supplied by --vina-bin
BOX = FIRST / "coarse_v1/config/vina_box.txt"
BOX_META = FIRST / "coarse_v1/config/common_box.json"
EXPORT = None  # supplied by --export-bin
SEED_PREP = 20260927
KEYS = ("standardized_isomeric_smiles", "rule_smiles", "smiles", "source_smiles")


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def read_rows(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_rows(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fields)
        w.writeheader()
        w.writerows(rows)


def get_input(row):
    compound_id = row.get("compound_id") or row.get("parent_id")
    if not compound_id:
        raise ValueError("missing compound_id")
    smiles = next((row[k] for k in KEYS if row.get(k)), "")
    if not smiles:
        raise ValueError("missing SMILES")
    return compound_id, smiles


def receptor_xyz(path):
    points = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_type = line.split()[-1]
        if atom_type in {"H", "HD", "HS"} or re.fullmatch(r"G[0-9]+", atom_type):
            continue
        points.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.asarray(points, dtype=float)


def pose_xyz(text):
    points = []
    for line in text.split("ENDMDL", 1)[0].splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_type = line.split()[-1]
        if atom_type in {"H", "HD", "HS"} or re.fullmatch(r"G[0-9]+", atom_type):
            continue
        points.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.asarray(points, dtype=float)


def prepare(item, ligdir):
    started = time.perf_counter()
    row = dict(item)
    cid, smiles = get_input(item)
    row["compound_id"] = cid
    row["preparation_started_utc"] = utc()
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("RDKit cannot parse input")
        if len(Chem.GetMolFrags(mol)) != 1:
            raise ValueError("multiple components")
        input_smi = Chem.MolToSmiles(mol, isomericSmiles=True)
        row["input_isomeric_smiles"] = input_smi
        row["input_inchikey"] = Chem.MolToInchiKey(mol)
        source_unknown = unspecified_potential_stereo(mol)
        row["source_unspecified_potential_stereo"] = json.dumps(source_unknown)
        if source_unknown:
            raise ValueError("input has unspecified potentially stable atom or double-bond stereochemistry")
        variants = sorted(set(protonate_smiles(input_smi, ph_min=7.4, ph_max=7.4,
                                               precision=0, max_variants=8)))
        row["rule_variant_count"] = len(variants)
        if len(variants) != 1:
            raise ValueError(f"rule pH7.4 returned {len(variants)} forms; no silent selection")
        state = Chem.MolFromSmiles(variants[0])
        if state is None:
            raise ValueError("invalid pH rule form")
        state_smi = Chem.MolToSmiles(state, isomericSmiles=True)
        row["rule_smiles"] = state_smi
        row["state_inchikey"] = Chem.MolToInchiKey(state)
        unspecified = unspecified_potential_stereo(state)
        exchangeable = exchangeable_n_indices(mol, state)
        row["rule_unspecified_stereocenters"] = json.dumps(unspecified)
        row["exchangeable_protonated_N_indices"] = json.dumps(exchangeable)
        row["chemical_identity_rule_version"] = RULE_VERSION
        if any(kind != "Atom_Tetrahedral" or idx not in exchangeable for kind,idx in unspecified):
            raise ValueError("rule pH7.4 form contains stable or unsupported unspecified stereocenter")
        row["exchangeable_N_geometry_approximation"] = bool(unspecified)
        row["heavy_atoms"] = state.GetNumHeavyAtoms()
        hm = Chem.AddHs(state)
        params = AllChem.ETKDGv3()
        params.randomSeed = SEED_PREP
        params.numThreads = 1
        params.pruneRmsThresh = 0.5
        confs = list(AllChem.EmbedMultipleConfs(hm, numConfs=10, params=params))
        row["conformers_generated"] = len(confs)
        if not confs or not AllChem.MMFFHasAllMoleculeParams(hm):
            raise ValueError("no conformer or MMFF parameters")
        en = AllChem.MMFFOptimizeMoleculeConfs(hm, numThreads=1, maxIters=1000,
                                               mmffVariant="MMFF94s")
        good = [(energy, conf) for conf, (status, energy) in zip(confs, en) if status == 0]
        row["conformers_converged"] = len(good)
        if not good:
            raise ValueError("no converged MMFF conformer")
        energy, conf = min(good)
        row["selected_mmff_energy_kcal_mol"] = float(energy)
        row["selected_conformer"] = int(conf)
        sdf = ligdir / f"{cid}.sdf"
        qt = ligdir / f"{cid}.pdbqt"
        with sdf.open("w", encoding="utf-8") as f:
            with Chem.SDWriter(f) as w:
                w.write(hm, confId=int(conf))
        with sdf.open("rb") as f:
            roundtrip = next(Chem.ForwardSDMolSupplier(f, removeHs=False))
        if roundtrip is None:
            raise ValueError("SDF readback failed")
        actual_mol = Chem.RemoveHs(roundtrip)
        after = Chem.MolToSmiles(actual_mol, isomericSmiles=True)
        row["sdf_graph_isomeric_smiles"] = after
        row["sdf_graph_identity_equal"] = after == state_smi
        row["sdf_graph_identity_stable_equal"] = stable_identity_smiles(actual_mol, exchangeable) == stable_identity_smiles(state, exchangeable)
        m3d = Chem.RemoveHs(roundtrip)
        Chem.AssignStereochemistryFrom3D(m3d, replaceExistingTags=True)
        row["sdf_3d_isomeric_smiles"] = Chem.MolToSmiles(m3d, isomericSmiles=True)
        row["sdf_3d_stereo_equal"] = row["sdf_3d_isomeric_smiles"] == state_smi
        row["sdf_3d_stable_identity_equal"] = stable_identity_smiles(m3d, exchangeable) == stable_identity_smiles(state, exchangeable)
        if not row["sdf_graph_identity_stable_equal"] or not row["sdf_3d_stable_identity_equal"]:
            raise ValueError("3D SDF stable identity/stereo mismatch; raw SDF retained")
        setups = MoleculePreparation().prepare(hm, conformer_id=int(conf))
        if len(setups) != 1:
            raise ValueError("not exactly one Meeko setup")
        pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            raise ValueError(err)
        qt.write_text(pdbqt, encoding="utf-8")
        row.update(preparation_status="prepared", sdf_file=str(sdf),
                   sdf_sha256=sha(sdf), ligand_pdbqt=str(qt), ligand_sha256=sha(qt))
    except Exception as e:
        row.update(preparation_status="blocked", preparation_reason=str(e),
                   preparation_traceback=traceback.format_exc())
    row["preparation_elapsed_seconds"] = round(time.perf_counter()-started, 3)
    return row


def dock(row, a, root, receptor, receptor_sha, box_sha, box_meta, tree):
    cid = row["compound_id"]
    run_id = f"{cid}__{a.receptor}__e{a.exhaustiveness}__s{a.seed}"
    pose = root / "poses" / f"{run_id}.pdbqt"
    sdf = root / "poses" / f"{run_id}.sdf"
    log = root / "logs" / f"{run_id}.log"
    export_log = root / "logs" / f"{run_id}_export.log"
    context = {"M_PUB_ALL": "mouse TMC1-CIB2-TMIE six-protein complex",
               "M_PUB_CORE": "mouse TMC1 core chain A only",
               "M_AF": "mouse TMC1 predicted structure, aligned",
               "H_AF": "human TMC1 predicted structure, aligned",
               "W_7USX_LOCAL_COMPLETE": "C. elegans experimental TMC1-CALM1-TMIE six-protein complex; remote incomplete residues omitted outside fixed box+8A"}
    result = dict(row, run_id=run_id, receptor_id=a.receptor,
                  protein_context=context[a.receptor],
                  lipids_present=False, conformation_evidence="7USX cryoEM contracted geometric conformation; functional closed/open state unverified",
                  binding_hypothesis="literature-informed middle-region pose",
                  functional_direction="unknown", disease_fit="unverified",
                  scoring_engine="AutoDock Vina 1.2.7", score_kind="Vina docking score, not MMGBSA",
                  seed=a.seed, exhaustiveness=a.exhaustiveness, num_modes=9,
                  energy_range=3.0, cpu_per_job=2,
                  receptor_file=str(receptor), receptor_sha256=receptor_sha,
                  box_file=str(a.box), box_sha256=box_sha,
                  box_id=box_meta["site_id"],
                  box_center_A=json.dumps(box_meta["center"]), box_size_A=json.dumps(box_meta["size"]),
                  pose_file=str(pose), log_file=str(log), export_log_file=str(export_log),
                  independent_qc_status="pending")
    command = [str(VINA), "--receptor", str(receptor), "--ligand", row["ligand_pdbqt"],
               "--config", str(a.box), "--scoring", "vina", "--exhaustiveness", str(a.exhaustiveness),
               "--num_modes", "9", "--energy_range", "3.0", "--cpu", "2",
               "--seed", str(a.seed), "--out", str(pose)]
    result["command"] = json.dumps(command)
    result["started_utc"] = utc()
    t = time.perf_counter()
    try:
        with log.open("w", encoding="utf-8") as f:
            proc = subprocess.run(command, stdout=f, stderr=subprocess.STDOUT,
                                  timeout=a.timeout)
        result["returncode"] = proc.returncode
        if proc.returncode:
            raise ValueError(f"Vina exit {proc.returncode}")
        text = pose.read_text(encoding="utf-8")
        scores = re.findall(r"REMARK VINA RESULT:\s*([-+0-9.]+)", text)
        if not scores:
            raise ValueError("pose has no Vina score")
        score = float(scores[0])
        if not math.isfinite(score):
            raise ValueError("nonfinite score")
        coords = pose_xyz(text)
        if len(coords) != int(row["heavy_atoms"]):
            raise ValueError("pose heavy atom count mismatch")
        distances = tree.query(coords, k=1)[0]
        center = np.asarray(box_meta["center"])
        half = np.asarray(box_meta["size"])/2
        margin = float(np.min(half-np.abs(coords-center)))
        result.update(status="completed", score_kcal_mol=score, pose_count=len(scores),
                      pose_sha256=sha(pose), pose_heavy_atoms=len(coords),
                      min_receptor_heavy_distance_A=float(distances.min()),
                      severe_heavy_atoms_below_1p5A=int((distances < 1.5).sum()),
                      min_box_margin_A=margin)
        with export_log.open("w", encoding="utf-8") as f:
            ex = subprocess.run([str(EXPORT), str(pose), "-s", str(sdf)],
                                stdout=f, stderr=subprocess.STDOUT, timeout=60)
        result["sdf_export_returncode"] = ex.returncode
        result["exported_sdf_file"] = str(sdf) if sdf.exists() else ""
        if sdf.exists():
            result["exported_sdf_sha256"] = sha(sdf)
        result["local_pre_qc_pass"] = bool(ex.returncode == 0 and sdf.exists() and
                                          distances.min() >= 1.5 and margin >= -0.1)
    except Exception as e:
        result.update(status="failed", failure_reason=str(e), local_pre_qc_pass=False,
                      traceback=traceback.format_exc())
    result["elapsed_seconds"] = round(time.perf_counter()-t, 3)
    result["finished_utc"] = utc()
    write_json(root / "logs" / f"{run_id}.json", result)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--batch", required=True)
    ap.add_argument("--scratch", type=Path, default=DEFAULT_SCRATCH)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--receptor", choices=["W_7USX_LOCAL_COMPLETE"], default="W_7USX_LOCAL_COMPLETE")
    ap.add_argument("--receptor-dir", required=True, type=Path)
    ap.add_argument("--vina-bin", required=True, type=Path)
    ap.add_argument("--export-bin", required=True, type=Path)
    ap.add_argument("--box", type=Path, default=BOX)
    ap.add_argument("--box-meta", type=Path, default=BOX_META)
    ap.add_argument("--seed", type=int, default=20260927)
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--timeout", type=int, default=1800)
    a = ap.parse_args()
    if not 1 <= a.workers <= 2:
        raise ValueError("maximum 2 workers /4 Vina CPU")
    global VINA, EXPORT
    VINA, EXPORT = a.vina_bin, a.export_bin
    if sha(VINA) != "e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5":
        raise ValueError("Vina binary SHA256 mismatch")
    if not EXPORT.is_file():
        raise FileNotFoundError("Meeko mk_export executable missing")
    if sha(a.box) != "17100e6a2f415d440b9f1aa74b3ecf2d53c681a0edd49d2b0573531512e4e149":
        raise ValueError("Vina box SHA256 mismatch")
    root = a.scratch / a.batch
    for name in ("ligands", "poses", "logs"):
        (root / name).mkdir(parents=True, exist_ok=True)
    receptor = a.receptor_dir / f"{a.receptor}.pdbqt"
    if sha(receptor) != "907ef675c3f92bbb3efbec7234cbbbf8a3fb1bd412f4191ce574e238ae9b6651":
        raise ValueError("Frozen 7USX receptor SHA256 mismatch")
    box_meta = json.loads(a.box_meta.read_text(encoding="utf-8"))
    input_rows = read_rows(a.input)
    if a.limit:
        input_rows = input_rows[:a.limit]
    seen = set()
    unique = []
    for r in input_rows:
        cid, _ = get_input(r)
        if cid not in seen:
            unique.append(r)
            seen.add(cid)
    plan = {"created_utc": utc(), "input_path": str(a.input), "input_sha256": sha(a.input),
            "input_rows": len(input_rows), "unique_ids": len(unique),
            "receptor": str(receptor), "receptor_sha256": sha(receptor),
            "box": box_meta, "box_file": str(a.box), "box_sha256": sha(a.box), "vina_sha256": sha(VINA),
            "rdkit": rdBase.rdkitVersion, "chemical_identity_rule_version": RULE_VERSION,
            "seed": a.seed, "conformer_seed": SEED_PREP,
            "exhaustiveness": a.exhaustiveness, "workers": a.workers,
            "protocol": "one pH7.4 rule form only; ETKDGv3 10 conformers, converged MMFF94s min; Vina first pose; no biological claim"}
    write_json(root / "plan.json", plan)
    prepared = []
    for i, r in enumerate(unique, 1):
        p = prepare(r, root / "ligands")
        prepared.append(p)
        write_json(root / "logs" / f"{p['compound_id']}_preparation.json", p)
        print("PREP", i, len(unique), p["compound_id"], p["preparation_status"], flush=True)
    write_rows(root / "preparation_ledger.csv", prepared)
    tree = cKDTree(receptor_xyz(receptor))
    jobs = []
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        fut = [ex.submit(dock, p, a, root, receptor, sha(receptor), sha(a.box), box_meta, tree)
               for p in prepared if p["preparation_status"] == "prepared"]
        for f in as_completed(fut):
            j = f.result()
            jobs.append(j)
            write_rows(root / "job_ledger.csv", sorted(jobs, key=lambda z: z["run_id"]))
            print("DOCK", len(jobs), len(fut), j["run_id"], j["status"],
                  j.get("score_kcal_mol", ""), j["elapsed_seconds"], flush=True)
    summary = {"finished_utc": utc(), "plan": plan, "prepared": len(prepared),
               "prepared_success": sum(p["preparation_status"] == "prepared" for p in prepared),
               "preparation_blocked": sum(p["preparation_status"] != "prepared" for p in prepared),
               "jobs": len(jobs), "completed": sum(j["status"] == "completed" for j in jobs),
               "local_pre_qc_pass": sum(j.get("local_pre_qc_pass", False) for j in jobs),
               "sum_job_seconds": sum(j["elapsed_seconds"] for j in jobs),
               "scratch_path": str(root)}
    write_json(root / "summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
