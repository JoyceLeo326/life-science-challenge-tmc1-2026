"""Mark retrospective docking-label overlap without modifying historical files."""
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import os

ROOT = Path(os.environ.get("TMC1_LIBRARY_WORK", str(Path(__file__).resolve().parent)))
FIRST = Path(os.environ.get("TMC1_R1_ROOT", str(Path(__file__).resolve().parent.parent / "provenance/historical")))
LIB = FIRST / "coarse_v1/inputs/library/selected_candidates.csv"
AUDIT = FIRST / "coarse_v1/analysis/job_audit_join.csv"
OUT = ROOT / "processed"


def rows(path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write(path, data, fields):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(data)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


historical = rows(LIB)
audited = rows(AUDIT)
hist_ids = {r["chembl_id"] for r in historical if r.get("chembl_id")}
qc_ids = {r["chembl_id"] for r in audited if r.get("final_qc_usable", "").lower() == "true"}
hist_keys = {r["rdkit_inchikey"].split("-")[0] for r in historical if r.get("rdkit_inchikey")}
qc_keys = {r["rdkit_inchikey"].split("-")[0] for r in historical if r.get("rdkit_inchikey") and r["chembl_id"] in qc_ids}
parents = rows(OUT / "all_unique_parents.csv")
eligible = rows(OUT / "filtered_for_screening.csv")


def flags(r):
    ids = set(r["source_ids"].split(";"))
    key = r["parent_key"]
    return {
        "historical_971_identity_overlap": bool(ids & hist_ids or key in hist_keys),
        "historical_935_qc_label_overlap": bool(ids & qc_ids or key in qc_keys),
        "historical_id_overlap": bool(ids & hist_ids),
        "historical_key_overlap": key in hist_keys,
    }


outrows = []
for r in parents:
    f = flags(r)
    outrows.append({"compound_id": r["compound_id"], "parent_key": r["parent_key"],
                    "source_ids": r["source_ids"], "processing_status": r["processing_status"],
                    **f})
write(OUT / "historical_overlap.csv", outrows, list(outrows[0]))
forward = [r for r in eligible if not flags(r)["historical_971_identity_overlap"]]
write(OUT / "forward_unlabeled_pool.csv", forward, list(eligible[0]))
counts = Counter()
for r in parents:
    f = flags(r)
    counts["all_historical_971_overlap"] += f["historical_971_identity_overlap"]
    counts["all_historical_935_qc_label_overlap"] += f["historical_935_qc_label_overlap"]
for r in eligible:
    f = flags(r)
    counts["eligible_historical_971_overlap"] += f["historical_971_identity_overlap"]
    counts["eligible_historical_935_qc_label_overlap"] += f["historical_935_qc_label_overlap"]
report = {"created_at_utc": datetime.now(timezone.utc).isoformat(), "historical_identity_rows": len(historical),
          "historical_qc_label_rows": len(qc_ids), "all_unique_parents": len(parents),
          "eligible_unique_parents": len(eligible), "forward_unlabeled_pool": len(forward),
          "overlap_counts": dict(counts), "matching_rule": "ChEMBL source ID OR nonstereo InChIKey first block",
          "historical_file_sha256": {str(LIB): sha(LIB), str(AUDIT): sha(AUDIT)}}
(OUT / "historical_overlap_summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(report, ensure_ascii=False), flush=True)
