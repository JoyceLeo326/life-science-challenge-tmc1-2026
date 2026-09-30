"""Build a read-only logical view of original 154 and resumed 124 jobs.

Call after both batches are quiescent or at a stable snapshot. It never edits
either producer batch. Job paths remain absolute to their source raw files.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from portable_config import required_path
from frozen_publication import frozen_matches

from rank_library import read_csv, write_csv, sha

HERE=required_path("TMC1_R2_ROOT") / "03_compute"
SCRATCH=required_path("TMC1_R2_ROOT") / "03_compute/raw_runs"
OLD=SCRATCH/"formal_forward_v2_initial278"
NEW=SCRATCH/"formal_forward_v2_resume124"
FROZEN=HERE/"ranking_formal_v2_16698"/"physical_union_docking_list.csv"
EXPECTED_FROZEN="3b53b86471eaed6f5bc49fbef381b811ec81153c69376fc4f20bf2c0161790b9"
OLD_LEDGER_SHA="75a95077a56ea1945b169b984d719b461620704c71b7e3d818771617173e14db"


def keyed(rows,kind):
    d={}
    for row in rows:
        cid=row["compound_id"]
        if cid in d:raise ValueError(f"duplicate {kind}: {cid}")
        d[cid]=row
    return d


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--require-complete",action="store_true")
    a=p.parse_args()
    if sha(FROZEN)!=EXPECTED_FROZEN or not frozen_matches(OLD/"job_ledger.csv",OLD_LEDGER_SHA):
        raise ValueError("original freeze or paused job ledger hash changed")
    frozen=read_csv(FROZEN)
    old_prep=keyed(read_csv(OLD/"preparation_ledger.csv"),"old prep")
    new_prep=keyed(read_csv(NEW/"preparation_ledger.csv"),"resume prep")
    old_jobs=keyed(read_csv(OLD/"job_ledger.csv"),"old jobs")
    new_jobs=keyed(read_csv(NEW/"job_ledger.csv"),"resume jobs")
    fids={r["compound_id"] for r in frozen}
    if len(frozen)!=278 or len(fids)!=278 or len(old_prep)!=278 or set(old_prep)!=fids:
        raise ValueError("original 278 prep/frozen identity mismatch")
    if len(old_jobs)!=154 or set(old_jobs)&set(new_jobs):
        raise ValueError("old and resumed jobs overlap or old count drift")
    remaining=fids-set(old_jobs)
    if len(remaining)!=124 or set(new_prep)!=remaining or not set(new_jobs)<=remaining:
        raise ValueError("resume prep/jobs differ from frozen remaining set")
    if a.require_complete and len(new_jobs)!=124:
        raise ValueError(f"resume incomplete: {len(new_jobs)}/124")
    all_prep={**{cid:old_prep[cid] for cid in old_jobs},**new_prep}
    all_jobs={**old_jobs,**new_jobs}
    if set(all_prep)!=fids:
        raise ValueError("merged prep not full frozen set")
    a.output.mkdir(parents=True,exist_ok=True)
    # Atomic replacement keeps readers from seeing a partial CSV.
    prep_tmp=a.output/"preparation_ledger.tmp.csv"
    job_tmp=a.output/"job_ledger.tmp.csv"
    write_csv(prep_tmp,[all_prep[r["compound_id"]] for r in frozen])
    write_csv(job_tmp,[all_jobs[r["compound_id"]] for r in frozen if r["compound_id"] in all_jobs])
    prep_tmp.replace(a.output/"preparation_ledger.csv")
    job_tmp.replace(a.output/"job_ledger.csv")
    manifest={"created_utc":datetime.now(timezone.utc).isoformat(),
              "status":"complete" if len(all_jobs)==278 else "partial_snapshot",
              "frozen_count":278,"prepared_count":len(all_prep),
              "old_completed_jobs":154,"resumed_completed_jobs":len(new_jobs),
              "merged_completed_jobs":len(all_jobs),
              "source_paths":{"old":str(OLD),"resumed":str(NEW)},
              "source_hashes":{"frozen":sha(FROZEN),
                               "old_jobs":sha(OLD/"job_ledger.csv"),
                               "resumed_jobs":sha(NEW/"job_ledger.csv"),
                               "old_prep":sha(OLD/"preparation_ledger.csv"),
                               "resumed_prep":sha(NEW/"preparation_ledger.csv")},
              "output_hashes":{"jobs":sha(a.output/"job_ledger.csv"),
                               "prep":sha(a.output/"preparation_ledger.csv")},
              "no_source_files_modified":True}
    plan={"created_utc":manifest["created_utc"],
          "input_path":str(FROZEN),"input_sha256":EXPECTED_FROZEN,
          "input_rows":278,"unique_ids":278,
          "logical_view":True,"source_batches":[OLD.name,NEW.name],
          "not_a_separate_physical_docking_run":True}
    (a.output/"plan.json").write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding="utf-8")
    if len(all_jobs)==278:
        summary={"finished_utc":manifest["created_utc"],"plan":plan,
                 "prepared":278,"jobs":278,
                 "completed":sum(j["status"]=="completed" for j in all_jobs.values()),
                 "logical_view":True,"source_batches":[str(OLD),str(NEW)]}
        (a.output/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    (a.output/"merge_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"status":manifest["status"],"merged_completed_jobs":len(all_jobs),
                      "resumed_completed_jobs":len(new_jobs),"output":str(a.output)},ensure_ascii=False))


if __name__=="__main__":
    main()
