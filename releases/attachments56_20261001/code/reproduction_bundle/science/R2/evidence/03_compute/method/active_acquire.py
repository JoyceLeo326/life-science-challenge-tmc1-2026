"""Frozen prospective second-stage ExtraTrees acquisition.

Consumes only the 60 AI-initial attempt outcomes, after independent QC. Other
arm docking scores are never an input to this script. Failed attempts retain
their budget cost but do not become regression labels.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sklearn.ensemble import ExtraTreesRegressor

from rank_library import FIRST, LABELS, OLD_LIB, SEED, FP_SIZE, feat, read_csv, write_csv, sha


def load_ai_attempts(path: Path, initial_ids: set[str]) -> tuple[list[dict], list[dict]]:
    rows = read_csv(path)
    if len(rows) != len(initial_ids):
        raise ValueError(f"expected {len(initial_ids)} AI initial rows, found {len(rows)}")
    ids = [r["compound_id"] for r in rows]
    if len(set(ids)) != len(ids) or set(ids) != initial_ids:
        raise ValueError("AI outcome IDs must equal the frozen AI initial IDs exactly")
    allowed = {"qc_pass", "preparation_failed", "docking_failed", "qc_failed"}
    for r in rows:
        if r["outcome"] not in allowed:
            raise ValueError(f"nonterminal or invalid outcome {r['compound_id']} {r['outcome']}")
        if r["outcome"] == "qc_pass":
            if not r.get("run_id") or not r.get("pose_sha256") or not r.get("independent_qc_evidence"):
                raise ValueError(f"QC-passing row lacks traceability {r['compound_id']}")
            score = float(r["score_kcal_mol"])
            if not math.isfinite(score):
                raise ValueError(f"nonfinite score {r['compound_id']}")
        elif r.get("score_kcal_mol"):
            raise ValueError(f"failed row must not contain training score {r['compound_id']}")
    return rows, [r for r in rows if r["outcome"] == "qc_pass"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--freeze", required=True, type=Path)
    ap.add_argument("--ai-attempts", required=True, type=Path,
                    help="Exactly 60 AI-initial terminal outcomes, independently QC-verified")
    ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    started = time.perf_counter()
    frozen = json.loads((a.freeze / "freeze.json").read_text(encoding="utf-8"))
    preflight = json.loads((a.freeze / "preflight_v2.json").read_text(encoding="utf-8"))
    for name, expected in preflight["files_sha256"].items():
        if sha(a.freeze / name) != expected:
            raise ValueError(f"frozen input changed: {name}")
    if sha(LABELS) != frozen["old_label_sha256"] or sha(OLD_LIB) != frozen["old_library_sha256"]:
        raise ValueError("historical training source hash changed")
    if frozen["AI_initial_queries"] != 60 or frozen["AI_active_remaining_queries"] != 60:
        raise ValueError("unexpected active learning budget")
    frozen_predictions = read_csv(a.freeze / "all_predictions.csv")
    frozen_assignments = read_csv(a.freeze / "arm_assignments.csv")
    initial = [r for r in frozen_assignments if r["arm"] == "AI_initial"]
    initial_ids = {r["compound_id"] for r in initial}
    initial_by_id = {r["compound_id"]: r for r in initial}
    if len(initial_ids) != 60:
        raise ValueError("frozen AI initial set invalid")
    outcome_manifest_path=a.ai_attempts.with_name("AI_initial_60_independent_outcomes_manifest.json")
    outcome_manifest=json.loads(outcome_manifest_path.read_text(encoding="utf-8"))
    if outcome_manifest["outcome_file_sha256"]!=sha(a.ai_attempts) or \
       outcome_manifest["ai_initial_assignments_sha256"]!=sha(a.freeze/"arm_assignments.csv"):
        raise ValueError("AI independent outcomes do not match their audit manifest and freeze")
    attempts, passed = load_ai_attempts(a.ai_attempts, initial_ids)
    passed.sort(key=lambda r: int(initial_by_id[r["compound_id"]]["selection_rank"]))
    warning_file=HERE/"formal_charge_warning_ids.csv"
    warning_ids={r["compound_id"] for r in read_csv(warning_file)}
    a.output.mkdir(parents=True, exist_ok=False)
    old_lib = {r["parent_id"]: r for r in read_csv(OLD_LIB)}
    training = [r for r in read_csv(LABELS)
                if r["status"] == "completed" and r.get("final_qc_usable") == "True"]
    if len(training) != frozen["model_training_qc_labels"]:
        raise ValueError("historical training label count drift")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=FP_SIZE)
    old_bits, old_desc, labels = [], [], []
    for r in training:
        _, _, bits, desc = feat(old_lib[r["parent_id"]]["rule_smiles"], gen)
        old_bits.append(bits)
        old_desc.append(desc)
        labels.append(float(r["score_kcal_mol"]))
    old_desc = np.asarray(old_desc, dtype=np.float32)
    med = np.median(old_desc, axis=0)
    scale = np.std(old_desc, axis=0)
    scale[scale == 0] = 1
    for r in passed:
        frozen_row = initial_by_id[r["compound_id"]]
        _, _, bits, desc = feat(frozen_row["rule_smiles"], gen)
        old_bits.append(bits)
        old_desc = np.vstack([old_desc, desc])
        labels.append(float(r["score_kcal_mol"]))
    matrix = np.hstack([np.asarray(old_bits, dtype=np.float32),
                        (old_desc-med)/scale])
    model = ExtraTreesRegressor(n_estimators=160, max_features=0.5,
                                min_samples_leaf=2, n_jobs=1,
                                random_state=SEED).fit(matrix, np.asarray(labels, dtype=np.float32))
    eligible = [r for r in frozen_predictions
                if r["docking_eligibility"] == "eligible" and r["compound_id"] not in initial_ids]
    features = []
    for r in eligible:
        _, _, bits, desc = feat(r["rule_smiles"], gen)
        features.append(np.concatenate([bits.astype(np.float32), (desc-med)/scale]))
    query_x = np.asarray(features, dtype=np.float32)
    tree_preds = np.asarray([tree.predict(query_x) for tree in model.estimators_], dtype=np.float32)
    means = tree_preds.mean(axis=0)
    sds = tree_preds.std(axis=0)
    scored = []
    for i, r in enumerate(eligible):
        scored.append({**r, "active_mean_kcal_mol": float(means[i]),
                       "active_tree_sd_kcal_mol": float(sds[i]),
                       "active_acquisition_kcal_mol": float(means[i] - 0.5*sds[i])})
    scored.sort(key=lambda r: (r["active_acquisition_kcal_mol"], r["compound_id"]))
    chosen = [{"arm": "AI_active", "selection_rank": 61+i, **r}
              for i, r in enumerate(scored[:60])]
    write_csv(a.output / "active_predictions_all_eligible.csv", scored)
    write_csv(a.output / "AI_active_60.csv", chosen)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_input": str(a.freeze),
        "frozen_file_sha256": {name: sha(a.freeze / name) for name in
                               ("all_predictions.csv", "arm_assignments.csv", "freeze.json")},
        "ai_attempts_file": str(a.ai_attempts),
        "ai_attempts_sha256": sha(a.ai_attempts),
        "ai_attempts_manifest_sha256": sha(outcome_manifest_path),
        "independent_qc_pass_training_count": len(passed),
        "charge_warning_ids_sha256": sha(warning_file),
        "qc_passing_charge_warning_training_ids": [r["compound_id"] for r in passed if r["compound_id"] in warning_ids],
        "ai_initial_failed_or_qc_failed_cost_count": len(attempts)-len(passed),
        "historical_training_count": len(training),
        "total_training_count": len(labels),
        "unqueried_eligible_count": len(eligible),
        "selected_count": len(chosen),
        "selection_rule": "mean ExtraTrees score minus 0.5 tree SD, ascending, compound ID ties",
        "labels_used": "935 historical QC scores plus independently QC-passing AI_initial scores only",
        "other_arms_actual_scores_read": False,
        "model_seed": SEED,
        "chemical_identity_rule_version": frozen["chemical_identity_rule_version"],
        "software": {"rdkit": rdBase.rdkitVersion,
                     "sklearn": __import__("sklearn").__version__},
        "output_sha256": {name: sha(a.output / name) for name in
                          ("active_predictions_all_eligible.csv", "AI_active_60.csv")},
        "wall_seconds": time.perf_counter()-started,
    }
    (a.output / "active_freeze.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
