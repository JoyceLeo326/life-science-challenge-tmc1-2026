"""Read-only 27+14 logical view of the frozen 41 new active physical jobs."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from portable_config import required_path
from frozen_publication import frozen_matches

from rank_library import read_csv, write_csv, sha

HERE=required_path("TMC1_R2_ROOT") / "03_compute"
SCRATCH=required_path("TMC1_R2_ROOT") / "03_compute/raw_runs"
OLD=SCRATCH/"formal_forward_v2_active_new"
NEW=SCRATCH/"formal_forward_v2_active_resume14"
FROZEN=HERE/"active_physical_v2"/"active_physical_new.csv"
RESUME=HERE/"method/active_settings_interrupt_checkpoint_20260929"/"active_resume_freeze"/"remaining_active_14.csv"
EXPECTED_FROZEN="24b4687a8f7421dc57a5edf5be3302388787f1e49c399a0c5afed7d7a35f31c7"
EXPECTED_RESUME="09ab9293ef2114cdf59041f6971a4e8318f14f264fbb00635cc8b1e2954855eb"
EXPECTED_OLD_JOBS="4c11c0ad5fffe2fc639503eb2424120846ce248e831db0f4d39db553996ca161"
EXPECTED_OLD_PREP="72c1a9b83d99ea507b1397fa8421d398c4b42817d802978deeeb4b60c4aa46ac"


def keyed(rows,kind):
    result={}
    for row in rows:
        cid=row["compound_id"]
        if cid in result:
            raise ValueError(f"duplicate {kind} ID {cid}")
        result[cid]=row
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if (OLD/"summary.json").exists() or not (NEW/"summary.json").exists():
        raise ValueError("source batches must be interrupted old and completed resumed")
    if sha(FROZEN)!=EXPECTED_FROZEN or sha(RESUME)!=EXPECTED_RESUME or \
       not frozen_matches(OLD/"job_ledger.csv",EXPECTED_OLD_JOBS) or not frozen_matches(OLD/"preparation_ledger.csv",EXPECTED_OLD_PREP):
        raise ValueError("frozen input or interrupted source changed")
    frozen=read_csv(FROZEN)
    remaining=read_csv(RESUME)
    old_prep=keyed(read_csv(OLD/"preparation_ledger.csv"),"old prep")
    new_prep=keyed(read_csv(NEW/"preparation_ledger.csv"),"new prep")
    old_jobs=keyed(read_csv(OLD/"job_ledger.csv"),"old job")
    new_jobs=keyed(read_csv(NEW/"job_ledger.csv"),"new job")
    ids=[r["compound_id"] for r in frozen]
    frozen_ids=set(ids)
    if len(ids)!=41 or len(frozen_ids)!=41 or set(old_prep)!=frozen_ids:
        raise ValueError("frozen/original preparation 41-set differs")
    if len(old_jobs)!=27 or len(new_jobs)!=14 or len(new_prep)!=14 or \
       set(old_jobs)&set(new_jobs) or set(old_jobs)|set(new_jobs)!=frozen_ids or \
       set(new_jobs)!=set(new_prep) or set(new_prep)!={r["compound_id"] for r in remaining}:
        raise ValueError("27+14 ID partition invalid")
    for cid,row in new_prep.items():
        old=old_prep[cid]
        for field in ("input_isomeric_smiles","input_inchikey","state_inchikey",
                      "rule_smiles","sdf_sha256","ligand_sha256","preparation_status"):
            if row.get(field)!=old.get(field):
                raise ValueError(f"resume preparation differs from original frozen state {cid} {field}")
    all_prep={**{cid:old_prep[cid] for cid in old_jobs},**new_prep}
    all_jobs={**old_jobs,**new_jobs}
    if any(j["status"]!="completed" for j in all_jobs.values()):
        raise ValueError("one or more physical jobs did not complete")
    args.output.mkdir(parents=True,exist_ok=False)
    write_csv(args.output/"preparation_ledger.csv",[all_prep[cid] for cid in ids])
    write_csv(args.output/"job_ledger.csv",[all_jobs[cid] for cid in ids])
    manifest={"created_utc":datetime.now(timezone.utc).isoformat(),
              "status":"complete","frozen_count":41,"old_completed_jobs":27,
              "resumed_completed_jobs":14,"merged_completed_jobs":41,
              "source_paths":{"old":str(OLD),"resumed":str(NEW)},
              "source_hashes":{"frozen":sha(FROZEN),"remaining":sha(RESUME),
                               "old_jobs":sha(OLD/"job_ledger.csv"),"resumed_jobs":sha(NEW/"job_ledger.csv"),
                               "old_prep":sha(OLD/"preparation_ledger.csv"),"resumed_prep":sha(NEW/"preparation_ledger.csv")},
              "output_hashes":{"jobs":sha(args.output/"job_ledger.csv"),"prep":sha(args.output/"preparation_ledger.csv")},
              "interruption_cause":"settings-switch exec session teardown; no scientific input changed",
              "no_source_files_modified":True,"no_duplicate_physical_queries":True}
    plan={"created_utc":manifest["created_utc"],"input_path":str(FROZEN),"input_sha256":EXPECTED_FROZEN,
          "input_rows":41,"unique_ids":41,"logical_view":True,
          "source_batches":[OLD.name,NEW.name],"not_a_separate_physical_docking_run":True}
    (args.output/"plan.json").write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding="utf-8")
    (args.output/"summary.json").write_text(json.dumps({"finished_utc":manifest["created_utc"],
          "plan":plan,"prepared":41,"jobs":41,"completed":41,"logical_view":True,
          "source_batches":[str(OLD),str(NEW)]},ensure_ascii=False,indent=2),encoding="utf-8")
    (args.output/"merge_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"merged_completed_jobs":41,"output":str(args.output),
                      "job_ledger_sha256":manifest["output_hashes"]["jobs"]},ensure_ascii=False))


if __name__=="__main__":
    main()
