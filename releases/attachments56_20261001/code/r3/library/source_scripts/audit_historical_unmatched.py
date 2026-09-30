"""Explain first-round identities absent from the current standardizable pool."""
import csv
from collections import defaultdict
from pathlib import Path
import os

ROOT = Path(os.environ.get("TMC1_LIBRARY_WORK", str(Path(__file__).resolve().parent)))
BASE = ROOT / "snapshots" / "v2_20000_new"
OLD = Path(os.environ.get("TMC1_R1_ROOT", str(Path(__file__).resolve().parent.parent / "provenance/historical"))) / "coarse_v1/inputs/library/selected_candidates.csv"


def rows(path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


history = rows(OLD)
source = rows(BASE / "all_source_records.csv")
parents = rows(BASE / "all_unique_parents.csv")
by_id = defaultdict(list)
for r in source:
    by_id[r["source_id"]].append(r)
present_ids = {sid for r in parents for sid in r["source_ids"].split(";")}
present_keys = {r["parent_key"] for r in parents}
missing = []
for r in history:
    sid = r["chembl_id"]
    key = r["rdkit_inchikey"].split("-")[0]
    if sid not in present_ids and key not in present_keys:
        missing.append({"chembl_id": sid, "old_parent_id": r["parent_id"],
                        "old_inchikey": r["rdkit_inchikey"],
                        "current_source_status": ";".join(x["processing_status"] + ":" + x["processing_reason"] for x in by_id[sid]) or "not_in_current_pages"})
target = ROOT / "reference_identity" / "historical_971_unmatched.csv"
with target.open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["chembl_id", "old_parent_id", "old_inchikey", "current_source_status"])
    w.writeheader()
    w.writerows(missing)
print(f"historical_unmatched={len(missing)} path={target}")
