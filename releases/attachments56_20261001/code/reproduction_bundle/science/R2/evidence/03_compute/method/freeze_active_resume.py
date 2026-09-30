"""Freeze the missing active physical IDs after the settings-switch interruption.

The original completed jobs are immutable; only absent IDs are rerun.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from portable_config import required_path
from frozen_publication import frozen_matches

from rank_library import read_csv, write_csv, sha

HERE=required_path("TMC1_R2_ROOT") / "03_compute"
INPUT=HERE/"active_physical_v2"/"active_physical_new.csv"
OLD=HERE/"raw_runs/formal_forward_v2_active_new"
EXPECTED_INPUT="24b4687a8f7421dc57a5edf5be3302388787f1e49c399a0c5afed7d7a35f31c7"
EXPECTED_OLD_JOBS="4c11c0ad5fffe2fc639503eb2424120846ce248e831db0f4d39db553996ca161"
EXPECTED_PREP="72c1a9b83d99ea507b1397fa8421d398c4b42817d802978deeeb4b60c4aa46ac"


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if sha(INPUT)!=EXPECTED_INPUT or not frozen_matches(OLD/"job_ledger.csv",EXPECTED_OLD_JOBS) or not frozen_matches(OLD/"preparation_ledger.csv",EXPECTED_PREP):
        raise ValueError("frozen active input or interrupted source ledger changed")
    if (OLD/"summary.json").exists():
        raise ValueError("interrupted source batch unexpectedly has a completed summary")
    source=read_csv(INPUT)
    prep=read_csv(OLD/"preparation_ledger.csv")
    jobs=read_csv(OLD/"job_ledger.csv")
    ids=[r["compound_id"] for r in source]
    completed={r["compound_id"] for r in jobs}
    if len(ids)!=41 or len(set(ids))!=41 or len(prep)!=41 or {r["compound_id"] for r in prep}!=set(ids):
        raise ValueError("active frozen/preparation 41-set mismatch")
    if len(jobs)!=27 or len(completed)!=27 or not completed<=set(ids) or any(r["status"]!="completed" for r in jobs):
        raise ValueError("interrupted completed job set mismatch")
    remaining=[r for r in source if r["compound_id"] not in completed]
    if len(remaining)!=14 or {r["compound_id"] for r in remaining}&completed:
        raise ValueError("remaining 14 partition invalid")
    args.output.mkdir(parents=True,exist_ok=False)
    out=args.output/"remaining_active_14.csv"
    write_csv(out,remaining)
    manifest={"created_utc":datetime.now(timezone.utc).isoformat(),
              "interruption_cause":"settings switch ended exec session; no scientific input change",
              "frozen_new_physical_count":41,"original_completed_count":27,
              "remaining_count":14,"all_original_completed_status":"completed",
              "input_path":str(INPUT),"input_sha256":sha(INPUT),
              "original_batch_path":str(OLD),
              "original_job_ledger_sha256":sha(OLD/"job_ledger.csv"),
              "original_preparation_ledger_sha256":sha(OLD/"preparation_ledger.csv"),
              "remaining_input_sha256":sha(out),
              "original_completed_ids":sorted(completed),
              "remaining_ids":[r["compound_id"] for r in remaining],
              "no_duplicate_physical_queries":True}
    (args.output/"resume_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({k:manifest[k] for k in ("original_completed_count","remaining_count","remaining_input_sha256","no_duplicate_physical_queries")},ensure_ascii=False))


if __name__=="__main__":
    main()
