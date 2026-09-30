"""Predict new-library Vina proxy scores and freeze a prospective three-arm query set.

Predictions are not docking or biological activity. Runs entirely from a hashed
library snapshot and first-round QC-passing docking labels. No new score is read.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from portable_config import required_path

HERE = Path(__file__).resolve().parent

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from dimorphite_dl import protonate_smiles
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import Ridge
from chemistry_rules import exchangeable_n_indices, unspecified_potential_stereo, RULE_VERSION

FIRST = required_path("TMC1_R1_ROOT")
LABELS = FIRST / "coarse_v1/analysis/job_audit_join.csv"
OLD_LIB = FIRST / "coarse_v1/inputs/library/selected_candidates.csv"
SEED = 20260929
FP_SIZE = 1024


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fields)
        w.writeheader()
        w.writerows(rows)


def feat(smiles, gen):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError("invalid smiles")
    fp = gen.GetFingerprint(m)
    arr = np.zeros(FP_SIZE, dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    d = np.asarray([Descriptors.MolWt(m), Descriptors.MolLogP(m),
                    Descriptors.TPSA(m), Descriptors.NumHDonors(m),
                    Descriptors.NumHAcceptors(m), Descriptors.NumRotatableBonds(m),
                    Descriptors.HeavyAtomCount(m), Chem.GetFormalCharge(m)], dtype=np.float32)
    return m, fp, arr, d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--per-arm", type=int, default=30)
    ap.add_argument("--ai-initial", type=int, default=0,
                    help="Initial AI queries before prospective active-learning update; 0 means full static arm")
    ap.add_argument("--exclude-ids", type=Path,
                    help="CSV with compound_id attempted in a separate pilot before this formal freeze")
    ap.add_argument("--expected-library-sha256", default="")
    a = ap.parse_args()
    start = time.perf_counter()
    a.output.mkdir(parents=True, exist_ok=True)
    library_sha = sha(a.library)
    if a.expected_library_sha256 and library_sha.lower() != a.expected_library_sha256.lower():
        raise ValueError("library hash changed")
    old_lib = {r["parent_id"]: r for r in read_csv(OLD_LIB)}
    jobs = read_csv(LABELS)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=FP_SIZE)
    known_parent_keys = {r["rdkit_inchikey"].split("-")[0] for r in old_lib.values() if r.get("rdkit_inchikey")}
    train_bits, train_desc, y, train_fps = [], [], [], []
    for r in jobs:
        if r["status"] != "completed" or r.get("final_qc_usable") != "True":
            continue
        _, fp, bits, d = feat(old_lib[r["parent_id"]]["rule_smiles"], gen)
        train_bits.append(bits)
        train_desc.append(d)
        train_fps.append(fp)
        y.append(float(r["score_kcal_mol"]))
    y = np.asarray(y, dtype=np.float32)
    td = np.asarray(train_desc, dtype=np.float32)
    med = np.median(td, axis=0)
    scale = np.std(td, axis=0)
    scale[scale == 0] = 1
    xtrain = np.hstack([np.asarray(train_bits, dtype=np.float32), (td-med)/scale])
    ridge = Ridge(alpha=100).fit((td-med)/scale, y)
    ai = ExtraTreesRegressor(n_estimators=160, max_features=0.5,
                             min_samples_leaf=2, n_jobs=1, random_state=SEED).fit(xtrain, y)
    source = read_csv(a.library)
    pilot_ids = {r["compound_id"] for r in read_csv(a.exclude_ids)} if a.exclude_ids else set()
    rows = []
    pred_bits, pred_desc, pred_fps, pred_indices = [], [], [], []
    counts = Counter()
    for src in source:
        r = dict(src)
        r["historical_label_overlap"] = False
        r["prediction_status"] = "pending"
        r["docking_eligibility"] = "pending"
        try:
            if src.get("processing_status") != "eligible":
                raise ValueError("source property filter excluded")
            mol, _, _, _ = feat(src["standardized_isomeric_smiles"], gen)
            source_unknown = unspecified_potential_stereo(mol)
            r["source_stereo_unspecified"] = bool(source_unknown)
            r["source_unspecified_potential_stereo"] = json.dumps(source_unknown)
            r["historical_label_overlap"] = src["parent_key"] in known_parent_keys
            variants = sorted(set(protonate_smiles(Chem.MolToSmiles(mol, isomericSmiles=True),
                                                   ph_min=7.4, ph_max=7.4,
                                                   precision=0, max_variants=8)))
            r["rule_variant_count"] = len(variants)
            if len(variants) != 1:
                raise ValueError("rule pH7.4 not unique")
            rule = Chem.MolFromSmiles(variants[0])
            if rule is None:
                raise ValueError("invalid rule form")
            r["rule_smiles"] = Chem.MolToSmiles(rule, isomericSmiles=True)
            r["rule_inchikey"] = Chem.MolToInchiKey(rule)
            rule_unknown = unspecified_potential_stereo(rule)
            exchangeable = exchangeable_n_indices(mol, rule)
            r["rule_stereo_unspecified"] = bool(rule_unknown)
            r["rule_unspecified_potential_stereo"] = json.dumps(rule_unknown)
            r["rule_exchangeable_N_indices"] = json.dumps(exchangeable)
            r["rule_stable_stereo_unspecified"] = any(kind != "Atom_Tetrahedral" or idx not in exchangeable for kind,idx in rule_unknown)
            r["exchangeable_N_geometry_approximation"] = bool(rule_unknown and not r["rule_stable_stereo_unspecified"])
            m, fp, bits, d = feat(r["rule_smiles"], gen)
            r["murcko_scaffold"] = MurckoScaffold.MurckoScaffoldSmiles(mol=m) or "ACYCLIC:"+Chem.MolToSmiles(m, isomericSmiles=False)
            r["prediction_status"] = "predicted"
            r["docking_eligibility"] = ("prepilot_already_attempted" if r["compound_id"] in pilot_ids else
                                        "previously_docked" if r["historical_label_overlap"] else
                                        "unspecified_source_stereo" if r["source_stereo_unspecified"] else
                                        "unspecified_rule_stable_stereo" if r["rule_stable_stereo_unspecified"] else "eligible")
            pred_indices.append(len(rows))
            pred_bits.append(bits)
            pred_desc.append(d)
            pred_fps.append(fp)
        except Exception as e:
            r["prediction_status"] = "failed"
            r["docking_eligibility"] = "blocked"
            r["prediction_failure"] = str(e)
        counts[r["docking_eligibility"]] += 1
        rows.append(r)
    if not pred_indices:
        raise RuntimeError("no molecules were predicable")
    pd = np.asarray(pred_desc, dtype=np.float32)
    x = np.hstack([np.asarray(pred_bits, dtype=np.float32), (pd-med)/scale])
    ridge_pred = ridge.predict((pd-med)/scale)
    tree_preds = np.asarray([t.predict(x) for t in ai.estimators_], dtype=np.float32)
    ai_pred = tree_preds.mean(axis=0)
    ai_sd = tree_preds.std(axis=0)
    nn_sim = np.asarray([max(DataStructs.BulkTanimotoSimilarity(fp, train_fps)) for fp in pred_fps])
    for local, ix in enumerate(pred_indices):
        rows[ix].update(descriptor_pred_kcal_mol=float(ridge_pred[local]),
                        ai_pred_kcal_mol=float(ai_pred[local]),
                        ai_tree_sd_kcal_mol=float(ai_sd[local]),
                        max_train_tanimoto=float(nn_sim[local]),
                        ood_flag_max_tanimoto_below_0p3=bool(nn_sim[local] < 0.3))
    write_csv(a.output / "all_predictions.csv", rows)
    eligible = [r for r in rows if r["docking_eligibility"] == "eligible"]
    if len(eligible) < 3*a.per_arm:
        raise RuntimeError("insufficient prospective pool")
    # All selectors are frozen now, before any new Vina result is opened.
    rng = np.random.default_rng(SEED)
    arm_lists = {
        "random": [eligible[i] for i in rng.permutation(len(eligible))[:a.per_arm]],
        "descriptor_ridge": sorted(eligible, key=lambda r: (r["descriptor_pred_kcal_mol"], r["compound_id"]))[:a.per_arm],
        "AI_initial": sorted(eligible, key=lambda r: (r["ai_pred_kcal_mol"], r["compound_id"]))[:a.ai_initial or a.per_arm],
    }
    arms = []
    union = {}
    for arm, group in arm_lists.items():
        for rank, r in enumerate(group, 1):
            arms.append({"arm": arm, "selection_rank": rank, **r})
            if r["compound_id"] not in union:
                union[r["compound_id"]] = dict(r, selected_arms=[])
            union[r["compound_id"]]["selected_arms"].append(arm)
    for r in union.values():
        r["selected_arms"] = ";".join(r["selected_arms"])
    write_csv(a.output / "arm_assignments.csv", arms)
    write_csv(a.output / "physical_union_docking_list.csv", list(union.values()))
    summary = {"created_utc": datetime.now(timezone.utc).isoformat(),
               "task": "prospective new-molecule proxy of old M_PUB_ALL e8 Vina score; no new docking used",
               "input": str(a.library), "input_sha256": library_sha, "old_label_sha256": sha(LABELS),
               "old_library_sha256": sha(OLD_LIB), "model_training_qc_labels": len(y),
               "source_rows": len(source), "prediction_status": dict(Counter(r["prediction_status"] for r in rows)),
               "docking_eligibility": dict(counts), "arm_budget_attempts": a.per_arm,
               "AI_initial_queries": len(arm_lists["AI_initial"]),
               "AI_active_remaining_queries": a.per_arm-len(arm_lists["AI_initial"]),
               "physical_union_initial": len(union),
               "arm_overlap_savings_initial": sum(len(g) for g in arm_lists.values())-len(union),
               "excluded_pilot_ids": len(pilot_ids),
               "excluded_pilot_file": str(a.exclude_ids) if a.exclude_ids else None,
               "excluded_pilot_sha256": sha(a.exclude_ids) if a.exclude_ids else None,
               "predeclared_favorable_score_threshold_kcal_mol": float(np.quantile(y, 0.1)),
               "threshold_source": "10th percentile of 935 historical QC-passing M_PUB_ALL e8 scores; docking-score endpoint only",
               "failure_policy": "Every selected attempted molecule consumes its arm budget even if preparation, execution, or QC fails; failures are not assigned a favorable score.",
               "primary_endpoint": "Among each arm's 120 attempted new parent molecules, number with finite QC-passing M_PUB_ALL e8 seed20260927 first-pose score <= frozen historical 10th percentile. Same attempt denominator even on failure; no biological activity interpretation.",
               "secondary_endpoint": "Within successful poses, predeclared middle-site heavy contact within 4A, reported separately from score endpoint; no post hoc replacement or reranking.",
               "active_learning_rule": "After AI initial batch only, refit same ExtraTrees on historical labels plus QC-passing AI-initial new scores; failed AI initial attempts retained as cost and not imputed. Select remaining by predicted mean minus 0.5 tree SD, deterministic compound ID ties, excluding already AI queried, no non-AI arm labels used.",
               "cost_policy": "Report per-arm sum actual preparation and docking job CPU-requested seconds plus wall-clock; physical overlapping jobs executed once, each arm charged individually.",
               "model_seed": SEED, "chemical_identity_rule_version": RULE_VERSION,
               "software": {"rdkit": rdBase.rdkitVersion, "sklearn": __import__("sklearn").__version__},
               "wall_seconds": time.perf_counter()-start}
    (a.output / "freeze.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
