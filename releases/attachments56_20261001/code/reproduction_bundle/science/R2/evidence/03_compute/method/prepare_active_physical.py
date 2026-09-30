"""Freeze the physical docking list for 60 AI-active attempts.

Any query already docked for random or Ridge is reused without looking at its
score. It still counts as an AI-active attempt in the final comparison.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from portable_config import required_path

from rank_library import read_csv, write_csv, sha


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--freeze",required=True,type=Path)
    p.add_argument("--active",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    preflight=json.loads((a.freeze/"preflight_v2.json").read_text(encoding="utf-8"))
    active_meta=json.loads((a.active/"active_freeze.json").read_text(encoding="utf-8"))
    if sha(a.freeze/"physical_union_docking_list.csv")!=preflight["files_sha256"]["physical_union_docking_list.csv"]:
        raise ValueError("initial physical freeze hash changed")
    if sha(a.active/"AI_active_60.csv")!=active_meta["output_sha256"]["AI_active_60.csv"]:
        raise ValueError("active selection hash changed")
    a.output.mkdir(parents=True,exist_ok=False)
    initial=read_csv(a.freeze/"physical_union_docking_list.csv")
    active=read_csv(a.active/"AI_active_60.csv")
    initial_ids={r["compound_id"] for r in initial}
    active_ids=[r["compound_id"] for r in active]
    if len(active)!=60 or len(set(active_ids))!=60:
        raise ValueError("active 60 invalid")
    new=[r for r in active if r["compound_id"] not in initial_ids]
    reused=[r for r in active if r["compound_id"] in initial_ids]
    if set(active_ids)!={r["compound_id"] for r in new+reused}:
        raise ValueError("partition mismatch")
    if new:
        write_csv(a.output/"active_physical_new.csv",new)
    if reused:
        write_csv(a.output/"active_reused_initial.csv",reused)
    manifest={"created_utc":datetime.now(timezone.utc).isoformat(),
              "selection_rule":"fixed active 60; prior physical overlap reused regardless of known score",
              "active_attempts":len(active),"new_physical_jobs":len(new),
              "reused_physical_jobs":len(reused),
              "final_physical_unique_attempts":len(initial)+len(new),
              "initial_union_sha256":sha(a.freeze/"physical_union_docking_list.csv"),
              "active_selection_sha256":sha(a.active/"AI_active_60.csv"),
              "new_physical_sha256":sha(a.output/"active_physical_new.csv") if new else None,
              "reused_sha256":sha(a.output/"active_reused_initial.csv") if reused else None,
              "no_new_scores_used_in_partition":True}
    (a.output/"active_physical_freeze.json").write_text(
        json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(manifest,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
