"""Deterministic, single-process matched-molecular-pair enumeration.

This is reaction-template/terminal-substituent enumeration, not a trained
generative model. Run with the first-round virtual environment Python.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from portable_config import required_path

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import AllChem, Descriptors, QED, rdMolDescriptors, rdFingerprintGenerator
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Contrib.SA_Score import sascorer

HERE = required_path("TMC1_RULE_OUTPUT")
BASE = required_path("TMC1_R1_ROOT") / "confirmation_v1"
INPUT = BASE / "inputs/library/selected_candidates.csv"
SCORES = BASE / "analysis/parent_setting_complete_denominators.csv"
PARENTS = ("CHEMBL2103870", "CHEMBL3894860")
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], keys: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def catalog(kind: FilterCatalogParams.FilterCatalogs) -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(kind)
    return FilterCatalog(params)


PAINS = catalog(FilterCatalogParams.FilterCatalogs.PAINS)
BRENK = catalog(FilterCatalogParams.FilterCatalogs.BRENK)


def alerts(mol: Chem.Mol, cat: FilterCatalog) -> list[str]:
    return sorted({hit.GetDescription() for hit in cat.GetMatches(mol)})


def properties(mol: Chem.Mol) -> dict:
    return {
        "mw": round(Descriptors.MolWt(mol), 3),
        "clogp": round(Descriptors.MolLogP(mol), 3),
        "tpsa_A2": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "formal_charge": Chem.GetFormalCharge(mol),
        "qed": round(QED.qed(mol), 4),
        "sa_score": round(sascorer.calculateScore(mol), 3),
        "pains": alerts(mol, PAINS),
        "brenk": alerts(mol, BRENK),
    }


def main() -> None:
    parents = {r["chembl_id"]: r for r in read_csv(INPUT) if r["chembl_id"] in PARENTS}
    assert set(parents) == set(PARENTS)
    score_rows = read_csv(SCORES)
    score_stats = {}
    for cid in PARENTS:
        values = sorted(float(r["median_score_kcal_mol"]) for r in score_rows
                        if r["parent_id"] == "LIB_" + cid
                        and r["all3_independent_qc"].lower() == "true")
        assert len(values) >= 5, (cid, values)
        score_stats[cid] = {"n_contexts": len(values), "median_kcal_mol": round((values[(len(values)-1)//2] + values[len(values)//2])/2, 3), "worst_kcal_mol": max(values)}

    attempted: list[dict] = []
    candidates: list[dict] = []

    def add(cid: str, label: str, method: str, rationale: str, smiles: str, route: str) -> None:
        parent_mol = Chem.MolFromSmiles(parents[cid]["source_smiles"])
        mol = Chem.MolFromSmiles(smiles)
        log = {"parent_id": "LIB_"+cid, "variant": label, "method": method, "reason": rationale, "route_hypothesis": route, "input_smiles": smiles}
        if mol is None:
            attempted.append({**log, "status": "invalid_smiles"})
            return
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        if canonical == Chem.MolToSmiles(parent_mol, isomericSmiles=True):
            attempted.append({**log, "status": "same_as_parent"})
            return
        sim = round(DataStructs.TanimotoSimilarity(FPGEN.GetFingerprint(mol), FPGEN.GetFingerprint(parent_mol)), 4)
        pp = properties(parent_mol)
        pr = properties(mol)
        reasons = []
        for key, low, high in (("mw", 150, 550), ("clogp", -1, 5), ("tpsa_A2", 30, 140), ("hbd", 0, 5), ("hba", 0, 10), ("rotatable_bonds", 0, 10), ("sa_score", 0, 5.5)):
            if not low <= pr[key] <= high:
                reasons.append(f"{key}_outside_{low}_{high}")
        if not 0.65 <= sim < 0.99:
            reasons.append("parent_ecfp4_similarity_outside_0.65_0.99")
        if set(pr["pains"]) - set(pp["pains"]):
            reasons.append("new_PAINS_alert")
        if set(pr["brenk"]) - set(pp["brenk"]):
            reasons.append("new_Brenk_alert")
        if pr["formal_charge"] != 0:
            reasons.append("non_neutral_source_structure")
        inchikey = Chem.MolToInchiKey(mol)
        compound_id = "DES_" + inchikey[:14]
        result = {
            "compound_id": compound_id,
            "parent_id": "LIB_"+cid,
            "parent_name": parents[cid]["name"],
            "parent_chembl_id": cid,
            "parent_source_url": parents[cid]["source_url"],
            "parent_source_smiles": parents[cid]["source_smiles"],
            "parent_inchikey": parents[cid]["source_standard_inchikey"],
            "design_method": method,
            "variant": label,
            "transformation_reason": rationale,
            "route_hypothesis": route,
            "original_candidate_smiles": smiles,
            "standardized_isomeric_smiles": canonical,
            "inchikey": inchikey,
            "parent_ecfp4_tanimoto": sim,
            "parent_median_vina_kcal_mol": score_stats[cid]["median_kcal_mol"],
            "parent_worst_vina_kcal_mol": score_stats[cid]["worst_kcal_mol"],
            "parent_qc_context_count": score_stats[cid]["n_contexts"],
            # Uncalibrated prior: 2 kcal/mol penalty at zero similarity.
            "docking_neighborhood_prior_kcal_mol": round(score_stats[cid]["median_kcal_mol"] + 2.0*(1-sim), 3),
            "mw": pr["mw"], "clogp": pr["clogp"], "tpsa_A2": pr["tpsa_A2"],
            "hbd": pr["hbd"], "hba": pr["hba"], "rotatable_bonds": pr["rotatable_bonds"],
            "heavy_atoms": pr["heavy_atoms"], "formal_charge": pr["formal_charge"],
            "qed": pr["qed"], "sa_score": pr["sa_score"],
            "pains_alerts": "|".join(pr["pains"]), "brenk_alerts": "|".join(pr["brenk"]),
            "parent_pains_alerts": "|".join(pp["pains"]), "parent_brenk_alerts": "|".join(pp["brenk"]),
            "constraint_status": "pass" if not reasons else "exclude",
            "constraint_reasons": "|".join(reasons),
            "public_exact_structure_status": "not_yet_searched",
            "evidence_class": "newly_enumerated_structure; no TMC1 bioactivity measured",
        }
        attempted.append({**log, "status": result["constraint_status"], "compound_id": compound_id, "reasons": result["constraint_reasons"]})
        candidates.append(result)

    # Lumacaftor has one terminal carboxylic acid. Amidation after activation
    # is a chemically recognizable operation, but synthetic feasibility and
    # binding geometry still need experimental confirmation.
    lum = parents["CHEMBL2103870"]["source_smiles"]
    acid = "C(=O)O"
    assert lum.count(acid) == 1
    amides = [
        ("primary_amide", "C(=O)N", "remove carboxylate charge while retaining carbonyl H-bonding", "NH3"),
        ("methylamide", "C(=O)NC", "probe small neutral amide substituent", "methylamine"),
        ("cyclopropylamide", "C(=O)NC1CC1", "probe compact hydrophobic occupancy at acid terminus", "cyclopropylamine"),
        ("hydroxyethylamide", "C(=O)NCCO", "recover terminal polarity after acid neutralization", "ethanolamine"),
        ("morpholine_amide", "C(=O)N1CCOCC1", "probe cyclic polar amide at terminus", "morpholine"),
    ]
    for label, replacement, rationale, amine in amides:
        add("CHEMBL2103870", label, "acid_to_amide_template_v1", rationale,
            lum.replace(acid, replacement),
            f"Parent benzoic acid -> activated acid -> coupling with {amine}; conceptual route, not executed")

    tri = parents["CHEMBL3894860"]["source_smiles"]
    old = "CN1CCN("
    assert tri.startswith(old) and tri.count(old) == 1
    n_substituents = [
        ("desmethyl", "N1CCN(", "compare terminal piperazine basicity and size", "N-dealkylation or assemble from N-unsubstituted piperazine precursor"),
        ("ethyl", "CCN1CCN(", "one-carbon homologation of N-methyl substituent", "alkylate N-unsubstituted piperazine precursor with ethyl electrophile"),
        ("hydroxyethyl", "OCCN1CCN(", "increase terminal aqueous polarity with small H-bond donor", "alkylate N-unsubstituted piperazine precursor with protected 2-hydroxyethyl electrophile"),
        ("acetyl", "CC(=O)N1CCN(", "neutralize distal piperazine nitrogen and add carbonyl acceptor", "acylate N-unsubstituted piperazine precursor with acetyl reagent"),
        ("cyclopropyl", "C1CC1N1CCN(", "probe compact lipophilic N-substitution", "N-arylation/reductive amination from N-unsubstituted piperazine precursor; route uncertain"),
    ]
    for label, replacement, rationale, route in n_substituents:
        add("CHEMBL3894860", label, "piperazine_N_substitution_template_v1", rationale,
            replacement + tri[len(old):], route + "; conceptual route, not executed")

    # Stable identity deduplication: stereochemical InChIKey, not conformer ID.
    seen = set()
    for r in candidates:
        if r["inchikey"] in seen:
            r["constraint_status"] = "exclude"
            r["constraint_reasons"] += "|duplicate_inchikey"
        seen.add(r["inchikey"])
    candidates.sort(key=lambda r: (r["parent_id"], r["variant"]))
    keys = list(candidates[0])
    write_csv(HERE / "design_constraints_all.csv", candidates, keys)
    write_csv(HERE / "transformation_log.csv", attempted,
              ["parent_id", "variant", "method", "reason", "route_hypothesis", "input_smiles", "status", "compound_id", "reasons"])

    passing = [r for r in candidates if r["constraint_status"] == "pass"]
    # Limit to four per parent, prioritizing polar or smaller analogs after
    # constraints. This list is fixed before any second-round docking.
    priority = {"primary_amide": 0, "methylamide": 1, "hydroxyethylamide": 2,
                "cyclopropylamide": 3, "morpholine_amide": 4,
                "desmethyl": 0, "hydroxyethyl": 1, "acetyl": 2, "ethyl": 3,
                "cyclopropyl": 4}
    passing.sort(key=lambda r: (r["parent_id"], priority[r["variant"]]))
    selected = []
    for cid in PARENTS:
        selected.extend([r for r in passing if r["parent_chembl_id"] == cid][:4])
    for idx, r in enumerate(selected, 1):
        r["fixed_list_order"] = idx
    fixed_keys = ["fixed_list_order"] + keys
    write_csv(HERE / "fixed_docking_list.csv", selected, fixed_keys)

    sdf_log = []
    # RDKit's Windows C++ file opener can reject non-ASCII paths; write the
    # MolBlock through Python's Unicode-capable file API instead.
    with (HERE / "preferred_candidates.sdf").open("w", encoding="utf-8", newline="\n") as handle:
        for r in selected:
            mol = Chem.AddHs(Chem.MolFromSmiles(r["standardized_isomeric_smiles"]))
            seed = 20260929 + r["fixed_list_order"]
            embed = AllChem.EmbedMolecule(mol, randomSeed=seed, maxAttempts=100, useRandomCoords=False)
            if embed == 0:
                try:
                    ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
                    ff.Minimize(maxIts=200)
                    energy = round(ff.CalcEnergy(), 3)
                except Exception:
                    energy = None
                mol.SetProp("_Name", r["compound_id"])
                handle.write(Chem.MolToMolBlock(mol))
                for k, v in {**r, "conformer_seed": seed, "MMFF_energy_kcal_mol": energy}.items():
                    handle.write(f">  <{k}>\n{v}\n\n")
                handle.write("$$$$\n")
                sdf_log.append({"compound_id": r["compound_id"], "seed": seed, "embed_status": "success", "MMFF_energy_kcal_mol": energy})
            else:
                sdf_log.append({"compound_id": r["compound_id"], "seed": seed, "embed_status": "failed", "MMFF_energy_kcal_mol": ""})
    write_csv(HERE / "sdf_generation_log.csv", sdf_log, list(sdf_log[0]))

    comparison = []
    for cid in PARENTS:
        parent = parents[cid]
        prop = properties(Chem.MolFromSmiles(parent["source_smiles"]))
        comparison.append({"compound_id": "LIB_"+cid, "role": "commercial_parent", "parent_id": "",
                           "name": parent["name"], "mw": prop["mw"], "clogp": prop["clogp"],
                           "tpsa_A2": prop["tpsa_A2"], "qed": prop["qed"], "sa_score": prop["sa_score"],
                           "parent_ecfp4_tanimoto": 1.0, "docking_neighborhood_prior_kcal_mol": score_stats[cid]["median_kcal_mol"],
                           "actual_docking_score_available": True})
    for r in candidates:
        comparison.append({"compound_id": r["compound_id"], "role": "enumerated_analogue", "parent_id": r["parent_id"],
                           "name": r["variant"], "mw": r["mw"], "clogp": r["clogp"],
                           "tpsa_A2": r["tpsa_A2"], "qed": r["qed"], "sa_score": r["sa_score"],
                           "parent_ecfp4_tanimoto": r["parent_ecfp4_tanimoto"],
                           "docking_neighborhood_prior_kcal_mol": r["docking_neighborhood_prior_kcal_mol"],
                           "actual_docking_score_available": False})
    write_csv(HERE / "parent_product_comparison_data.csv", comparison, list(comparison[0]))

    manifest = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(), "method": "deterministic reaction/terminal-substituent enumeration",
        "is_trained_generative_model": False, "rdkit_version": rdBase.rdkitVersion,
        "input_paths_and_sha256": {str(p): sha256(p) for p in (INPUT, SCORES)},
        "parent_ids": list(PARENTS), "attempted_transformations": len(attempted),
        "valid_unique_structures": len(candidates), "constraint_pass": len(passing),
        "fixed_for_03_line": len(selected), "sdf_success": sum(x["embed_status"] == "success" for x in sdf_log),
        "proxy_definition": "median parent Vina across independently QC-passing receptor/representation settings + 2.0*(1-Morgan radius2/2048 Tanimoto); uncalibrated neighborhood prior, NOT a measured or predicted affinity",
        "design_status": "in silico enumerated only; no binding, function, synthesis or availability demonstrated",
    }
    (HERE / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    os.environ["OMP_NUM_THREADS"] = "1"
    main()
