"""Locked-budget prospective docking-score comparison after independent QC.

Failure consumes a query. Docking and contact are separate endpoints. Every
reported favorable score requires a matching independent QC pass record.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

from rank_library import OLD_LIB, LABELS, read_csv, write_csv, sha

HERE=Path(__file__).resolve().parent
CHARGE_WARNINGS=HERE/"formal_charge_warning_ids.csv"
STANDBY=HERE/"standby_intervals_20260929.csv"
INTERRUPTED_ACTIVE=HERE/"active_settings_interrupt_checkpoint_20260929"/"interrupted_inflight_attempts.csv"


def standby_intervals():
    return [(datetime.fromisoformat(r["start_local"]),
             datetime.fromisoformat(r["end_local"])) for r in read_csv(STANDBY)]


def standby_overlap_seconds(start: str, finish: str) -> float:
    if not start or not finish:
        return 0.0
    a=datetime.fromisoformat(start)
    b=datetime.fromisoformat(finish)
    if a.tzinfo is None or b.tzinfo is None or b<a:
        raise ValueError("job timestamps must be ordered and timezone-aware")
    return round(sum(max(0.0,(min(b,end)-max(a,begin)).total_seconds())
                     for begin,end in standby_intervals()),3)


def awake_job_seconds(job: dict) -> float:
    if not job or not job.get("elapsed_seconds"):
        return 0.0
    raw=float(job["elapsed_seconds"])
    corrected=raw-standby_overlap_seconds(job.get("started_utc",""),
                                           job.get("finished_utc",""))
    if corrected < -2.0:
        raise ValueError("standby overlap exceeds recorded docking duration")
    return round(max(0.0,corrected),3)


def to_bool(x: str) -> bool:
    if x not in {"True", "False"}:
        raise ValueError(f"expected explicit True or False, got {x!r}")
    return x == "True"


def by_id(rows, kind):
    result={}
    for r in rows:
        cid=r["compound_id"]
        if cid in result:
            raise ValueError(f"duplicate {kind} compound {cid}")
        result[cid]=r
    return result


def finite_or_none(value):
    if value in (None, ""):
        return None
    result=float(value)
    return result if math.isfinite(result) else None


def batch_wall_seconds(root: Path, old_override: Path | None = None, resumed_override: Path | None = None) -> float:
    merged=root/"merge_manifest.json"
    if merged.exists():
        meta=json.loads(merged.read_text(encoding="utf-8"))
        old=old_override or Path(meta["source_paths"]["old"])
        new=resumed_override or Path(meta["source_paths"]["resumed"])
        old_plan=json.loads((old/"plan.json").read_text(encoding="utf-8"))
        old_jobs=read_csv(old/"job_ledger.csv")
        old_start=datetime.fromisoformat(old_plan["created_utc"])
        old_stop=max(datetime.fromisoformat(j["finished_utc"]) for j in old_jobs)
        return round((old_stop-old_start).total_seconds()+batch_wall_seconds(new),3)
    summary=json.loads((root/"summary.json").read_text(encoding="utf-8"))
    start=datetime.fromisoformat(summary["plan"]["created_utc"])
    finish=datetime.fromisoformat(summary["finished_utc"])
    return round((finish-start).total_seconds(),3)


def batch_awake_wall_seconds(root: Path, old_override: Path | None = None, resumed_override: Path | None = None) -> float:
    merged=root/"merge_manifest.json"
    if merged.exists():
        meta=json.loads(merged.read_text(encoding="utf-8"))
        old=old_override or Path(meta["source_paths"]["old"])
        new=resumed_override or Path(meta["source_paths"]["resumed"])
        old_plan=json.loads((old/"plan.json").read_text(encoding="utf-8"))
        old_jobs=read_csv(old/"job_ledger.csv")
        old_start=old_plan["created_utc"]
        old_stop=max(old_jobs,key=lambda j:j["finished_utc"])["finished_utc"]
        return round((datetime.fromisoformat(old_stop)-datetime.fromisoformat(old_start)).total_seconds()
                     -standby_overlap_seconds(old_start,old_stop)+batch_awake_wall_seconds(new),3)
    summary=json.loads((root/"summary.json").read_text(encoding="utf-8"))
    start=summary["plan"]["created_utc"]
    finish=summary["finished_utc"]
    return round(batch_wall_seconds(root)-standby_overlap_seconds(start,finish),3)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--freeze",required=True,type=Path)
    p.add_argument("--active",required=True,type=Path)
    p.add_argument("--initial-old-batch",required=True,type=Path)
    p.add_argument("--initial-resumed-batch",required=True,type=Path)
    p.add_argument("--active-old-batch",required=True,type=Path)
    p.add_argument("--active-resumed-batch",required=True,type=Path)
    p.add_argument("--initial-batch",required=True,type=Path)
    p.add_argument("--active-batch",type=Path,
                   help="Omit only if all active selections overlap initial physical queries")
    p.add_argument("--initial-qc",required=True,type=Path)
    p.add_argument("--active-qc",type=Path)
    p.add_argument("--initial-contacts",required=True,type=Path)
    p.add_argument("--active-contacts",type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    for batch in (a.initial_old_batch, a.initial_resumed_batch, a.active_old_batch, a.active_resumed_batch):
        if not (batch/"plan.json").is_file() or not (batch/"job_ledger.csv").is_file():
            raise FileNotFoundError(f"Physical batch plan/job ledger missing: {batch}")
    meta=json.loads((a.freeze/"freeze.json").read_text(encoding="utf-8"))
    preflight=json.loads((a.freeze/"preflight_v2.json").read_text(encoding="utf-8"))
    active_meta=json.loads((a.active/"active_freeze.json").read_text(encoding="utf-8"))
    for name,expected in preflight["files_sha256"].items():
        if sha(a.freeze/name)!=expected:
            raise ValueError(f"formal freeze file changed: {name}")
    if sha(a.active/"AI_active_60.csv")!=active_meta["output_sha256"]["AI_active_60.csv"]:
        raise ValueError("AI active selection changed")
    if sha(a.freeze/"arm_assignments.csv")!=active_meta["frozen_file_sha256"]["arm_assignments.csv"]:
        raise ValueError("AI active used a different initial freeze")
    initial_plan=json.loads((a.initial_batch/"plan.json").read_text(encoding="utf-8"))
    if initial_plan["input_sha256"]!=preflight["files_sha256"]["physical_union_docking_list.csv"]:
        raise ValueError("initial Vina batch input differs from formal freeze")
    threshold=float(meta["predeclared_favorable_score_threshold_kcal_mol"])
    warning_ids={r["compound_id"] for r in read_csv(CHARGE_WARNINGS)}
    if warning_ids!={"LIB_PLFCIUGPHZFBAC","LIB_DIMKSBQDWMXKAH"}:
        raise ValueError("charge warning list drift")
    arms=read_csv(a.freeze/"arm_assignments.csv")+read_csv(a.active/"AI_active_60.csv")
    if Counter(r["arm"] for r in arms)!={"random":120,"descriptor_ridge":120,"AI_initial":60,"AI_active":60}:
        raise ValueError("arm budgets differ from frozen plan")
    initial_ids={r["compound_id"] for r in read_csv(a.freeze/"physical_union_docking_list.csv")}
    active_ids={r["compound_id"] for r in arms if r["arm"]=="AI_active"}
    if len(active_ids)!=60 or active_ids & {r["compound_id"] for r in arms if r["arm"]=="AI_initial"}:
        raise ValueError("AI active repeats an AI initial attempt")
    a.output.mkdir(parents=True,exist_ok=False)
    initial_prep=by_id(read_csv(a.initial_batch/"preparation_ledger.csv"),"initial prep")
    initial_jobs=by_id(read_csv(a.initial_batch/"job_ledger.csv"),"initial job")
    if set(initial_prep)!=initial_ids:
        raise ValueError("initial batch attempts differ from frozen physical union")
    new_active_ids=active_ids-initial_ids
    interrupted_active=read_csv(INTERRUPTED_ACTIVE) if (a.active_batch and (a.active_batch/"merge_manifest.json").exists()) else []
    if interrupted_active and (len(interrupted_active)!=4 or
                               not {r["compound_id"] for r in interrupted_active}<=new_active_ids or
                               any(r["terminal_job_record"]!="False" or r["resource_time"]!="unknown"
                                   for r in interrupted_active)):
        raise ValueError("settings-switch interrupted attempt evidence mismatch")
    if new_active_ids:
        if not a.active_batch or not a.active_qc or not a.active_contacts:
            raise ValueError("active new queries require batch, QC and contacts")
        active_prep=by_id(read_csv(a.active_batch/"preparation_ledger.csv"),"active prep")
        active_jobs=by_id(read_csv(a.active_batch/"job_ledger.csv"),"active job")
        if set(active_prep)!=new_active_ids:
            raise ValueError("active batch attempts differ from frozen new active union")
    else:
        active_prep={}; active_jobs={}
    preps={**initial_prep,**active_prep}
    jobs={**initial_jobs,**active_jobs}
    qc=by_id(read_csv(a.initial_qc)+(read_csv(a.active_qc) if a.active_qc else []),"QC")
    contacts=by_id(read_csv(a.initial_contacts)+(read_csv(a.active_contacts) if a.active_contacts else []),"contacts")
    old_scaffolds=set()
    old_by_id={r["parent_id"]:r for r in read_csv(OLD_LIB)}
    historical_training=[r for r in read_csv(LABELS)
                         if r["status"]=="completed" and r.get("final_qc_usable")=="True"]
    if len(historical_training)!=935:
        raise ValueError("historical QC training set drift")
    for hist in historical_training:
        old=old_by_id[hist["parent_id"]]
        mol=Chem.MolFromSmiles(old["rule_smiles"])
        old_scaffolds.add(MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or
                          "ACYCLIC:"+Chem.MolToSmiles(mol,isomericSmiles=False))
    out=[]
    for r in arms:
        cid=r["compound_id"]
        prep=preps[cid]
        if prep.get("rule_smiles")!=r["rule_smiles"]:
            raise ValueError(f"preparation rule form differs from frozen selection {cid}")
        job=jobs.get(cid)
        independent=qc.get(cid)
        c=contacts.get(cid)
        prep_ok=prep["preparation_status"]=="prepared"
        if prep_ok and job is None:
            raise ValueError(f"prepared query has no terminal docking job {cid}")
        if job and job["status"]=="completed" and independent is None:
            raise ValueError(f"completed query lacks independent QC {cid}")
        qcp=to_bool(independent["independent_qc_pass"]) if independent else False
        if qcp and (not job or job["status"]!="completed" or
                    independent["run_id"]!=job["run_id"] or
                    independent["pose_sha256"]!=job["pose_sha256"]):
            raise ValueError(f"QC pass mismatch {cid}")
        if qcp and (c is None or c["run_id"]!=job["run_id"] or
                    c["pose_sha256"]!=job["pose_sha256"]):
            raise ValueError(f"QC-passing pose lacks matching region contacts {cid}")
        score=finite_or_none(job.get("score_kcal_mol")) if qcp else None
        if qcp and score is None:
            raise ValueError(f"QC passed but score not finite {cid}")
        if qcp and (finite_or_none(independent.get("score_kcal_mol")) is None or
                    abs(score-float(independent["score_kcal_mol"]))>1e-6):
            raise ValueError(f"independent QC score differs from producer ledger {cid}")
        contact=to_bool(c["any_site_contact_4A"]) if qcp else False
        if qcp and independent.get("middle_contact_4A") in {"True","False"}:
            if contact!=to_bool(independent["middle_contact_4A"]):
                raise ValueError(f"independent and local middle contacts disagree {cid}")
        if qcp and independent.get("closest_middle_A"):
            if abs(float(c["min_site_distance_A"])-float(independent["closest_middle_A"]))>0.02:
                raise ValueError(f"independent and local middle distance disagree {cid}")
        mol=Chem.MolFromSmiles(r["rule_smiles"])
        charge=Chem.GetFormalCharge(mol)
        actual_failure=("qc_pass" if qcp else "preparation_failed" if not prep_ok else
                        "docking_failed" if job["status"]!="completed" else "independent_qc_failed")
        out.append({
            "arm":r["arm"],"selection_rank":r["selection_rank"],"compound_id":cid,
            "attempt_outcome":actual_failure,
            "physical_batch":a.initial_batch.name if cid in initial_ids else a.active_batch.name,
            "run_id":job["run_id"] if job else "",
            "independent_qc_pass":qcp,
            "score_kcal_mol":score if score is not None else "",
            "favorable_score_predeclared":bool(qcp and score<=threshold),
            "middle_contact_4A":bool(qcp and contact),
            "favorable_and_middle_contact":bool(qcp and score<=threshold and contact),
            "min_middle_distance_A":c.get("min_site_distance_A","") if c else "",
            "min_box_margin_A":job.get("min_box_margin_A","") if job else "",
            "max_train_tanimoto":r.get("max_train_tanimoto",""),
            "ood_tanimoto_below_0p3":r.get("ood_flag_max_tanimoto_below_0p3",""),
            "murcko_scaffold":r.get("murcko_scaffold",""),
            "historical_scaffold_seen":r.get("murcko_scaffold","") in old_scaffolds,
            "rule_formal_charge":charge,"rule_abs_charge_ge2":abs(charge)>=2,
            "Gasteiger_charge_warning":cid in warning_ids,
            "rule_mw":round(Descriptors.MolWt(mol),3),
            "preparation_elapsed_seconds":prep.get("preparation_elapsed_seconds",""),
            "docking_wall_seconds":job.get("elapsed_seconds","") if job else "",
            "docking_standby_overlap_seconds":standby_overlap_seconds(job.get("started_utc",""),job.get("finished_utc","")) if job else "",
            "docking_awake_wall_seconds":awake_job_seconds(job) if job else "",
            "docking_requested_cpu_seconds":2*awake_job_seconds(job) if job else "",
            "independent_qc_evidence":independent.get("independent_qc_evidence","") if independent else "",
        })
    write_csv(a.output/"forward_attempt_level.csv",out)
    grouped=defaultdict(list)
    for r in out:
        arm="AI" if r["arm"].startswith("AI_") else r["arm"]
        grouped[arm].append(r)
    if {k:len(v) for k,v in grouped.items()}!={"random":120,"descriptor_ridge":120,"AI":120}:
        raise ValueError("final arm denominators differ")
    summary=[]
    for arm,rows in grouped.items():
        rec={"arm":arm,"attempted":len(rows),"qc_pass":sum(r["independent_qc_pass"] for r in rows),
             "preparation_failed":sum(r["attempt_outcome"]=="preparation_failed" for r in rows),
             "docking_failed":sum(r["attempt_outcome"]=="docking_failed" for r in rows),
             "independent_qc_failed":sum(r["attempt_outcome"]=="independent_qc_failed" for r in rows),
             "favorable_score":sum(r["favorable_score_predeclared"] for r in rows),
             "favorable_score_per_attempt":sum(r["favorable_score_predeclared"] for r in rows)/len(rows),
             "middle_contact_4A":sum(r["middle_contact_4A"] for r in rows),
             "favorable_and_middle_contact":sum(r["favorable_and_middle_contact"] for r in rows),
             "rule_abs_charge_ge2":sum(r["rule_abs_charge_ge2"] for r in rows),
             "historical_scaffold_seen":sum(r["historical_scaffold_seen"] for r in rows),
             "ood_tanimoto_below_0p3":sum(r["ood_tanimoto_below_0p3"]=="True" for r in rows),
             "sum_preparation_wall_seconds":round(sum(float(r["preparation_elapsed_seconds"] or 0) for r in rows),3),
             "sum_docking_wall_seconds":round(sum(float(r["docking_wall_seconds"] or 0) for r in rows),3),
             "sum_docking_awake_wall_seconds":round(sum(float(r["docking_awake_wall_seconds"] or 0) for r in rows),3),
             "sum_docking_requested_cpu_seconds":round(sum(float(r["docking_requested_cpu_seconds"] or 0) for r in rows),3)}
        summary.append(rec)
    write_csv(a.output/"forward_arm_summary.csv",summary)
    sensitivity=[]
    for arm,rows in grouped.items():
        flagged=[r for r in rows if r["Gasteiger_charge_warning"]]
        sensitivity.append({
            "arm":arm,"attempted_fixed":len(rows),
            "charge_warning_attempts":len(flagged),
            "primary_qc_pass":sum(r["independent_qc_pass"] for r in rows),
            "primary_favorable":sum(r["favorable_score_predeclared"] for r in rows),
            "qc_pass_if_warnings_failed":sum(r["independent_qc_pass"] for r in rows if not r["Gasteiger_charge_warning"]),
            "favorable_if_warnings_failed":sum(r["favorable_score_predeclared"] for r in rows if not r["Gasteiger_charge_warning"]),
            "interpretation":"conditional on the actually selected arm; not a rerun of active acquisition"})
    write_csv(a.output/"forward_charge_sensitivity.csv",sensitivity)
    strata=[]
    for arm,rows in grouped.items():
        criteria={"known_scaffold":lambda r:r["historical_scaffold_seen"],
                  "new_scaffold":lambda r:not r["historical_scaffold_seen"],
                  "tanimoto_below_0p3":lambda r:r["ood_tanimoto_below_0p3"]=="True",
                  "tanimoto_at_least_0p3":lambda r:r["ood_tanimoto_below_0p3"]!="True",
                  "abs_charge_ge2":lambda r:r["rule_abs_charge_ge2"],
                  "abs_charge_le1":lambda r:not r["rule_abs_charge_ge2"],
                  "mw_below_350":lambda r:r["rule_mw"]<350,
                  "mw_at_least_350":lambda r:r["rule_mw"]>=350}
        for label,predicate in criteria.items():
            selected=[r for r in rows if predicate(r)]
            strata.append({"arm":arm,"stratum":label,"attempted":len(selected),
                           "qc_pass":sum(r["independent_qc_pass"] for r in selected),
                           "favorable_score":sum(r["favorable_score_predeclared"] for r in selected),
                           "middle_contact_4A":sum(r["middle_contact_4A"] for r in selected),
                           "favorable_and_middle_contact":sum(r["favorable_and_middle_contact"] for r in selected)})
    write_csv(a.output/"forward_strata.csv",strata)
    audit={"created_utc":datetime.now(timezone.utc).isoformat(),
           "threshold_kcal_mol":threshold,
           "attempt_denominator_per_arm":120,
           "physical_initial_queries":len(initial_ids),
           "physical_new_active_queries":len(new_active_ids),
           "physical_total_unique_queries":len(initial_ids|new_active_ids),
           "arm_overlap_savings":360-len(initial_ids|new_active_ids),
           "physical_sum_preparation_wall_seconds":round(sum(float(r.get("preparation_elapsed_seconds") or 0) for r in preps.values()),3),
           "physical_sum_docking_wall_seconds":round(sum(float(r.get("elapsed_seconds") or 0) for r in jobs.values()),3),
           "physical_sum_docking_awake_wall_seconds":round(sum(awake_job_seconds(r) for r in jobs.values()),3),
           "physical_sum_docking_requested_cpu_seconds":round(2*sum(awake_job_seconds(r) for r in jobs.values()),3),
           "initial_batch_wall_seconds":batch_wall_seconds(a.initial_batch, a.initial_old_batch, a.initial_resumed_batch),
           "initial_batch_awake_wall_seconds":batch_awake_wall_seconds(a.initial_batch, a.initial_old_batch, a.initial_resumed_batch),
           "initial_batch_wall_excludes_user_pause":(a.initial_batch/"merge_manifest.json").exists(),
           "active_batch_wall_seconds":batch_wall_seconds(a.active_batch, a.active_old_batch, a.active_resumed_batch) if new_active_ids else 0,
           "active_batch_awake_wall_seconds":batch_awake_wall_seconds(a.active_batch, a.active_old_batch, a.active_resumed_batch) if new_active_ids else 0,
           "active_batch_wall_excludes_settings_switch_gap":bool(new_active_ids and (a.active_batch/"merge_manifest.json").exists()),
           "interrupted_inflight_active_attempts":len(interrupted_active),
           "interrupted_inflight_active_attempt_resource_seconds":None if interrupted_active else 0,
           "interrupted_inflight_manifest_sha256":sha(INTERRUPTED_ACTIVE) if interrupted_active else None,
           "standby_intervals_sha256":sha(STANDBY),
           "timing_note":"observed wall contains system standby; awake sums subtract standby overlap; merged batches exclude user pause or settings-switch interruption; interrupted unfinished docking work has unknown resource cost and is excluded from terminal job sums; 2 CPU times awake wall is slot capacity estimate, not measured CPU time",
           "shared_static_model_and_prediction_wall_seconds":meta.get("wall_seconds"),
           "active_refit_and_acquisition_wall_seconds":active_meta.get("wall_seconds"),
           "source_hashes":{str(path):sha(path) for path in
                            [a.freeze/"arm_assignments.csv",a.active/"AI_active_60.csv",
                             a.initial_batch/"preparation_ledger.csv",a.initial_batch/"job_ledger.csv",
                             a.initial_qc,a.initial_contacts]+
                            ([a.active_batch/"preparation_ledger.csv",a.active_batch/"job_ledger.csv",
                              a.active_qc,a.active_contacts] if new_active_ids else [])},
           "charge_warning_list_sha256":sha(CHARGE_WARNINGS),
           "interpretation":"Vina docking-score and pose geometry only; not activity, efficacy, or channel direction",
           "output_hashes":{name:sha(a.output/name) for name in
                            ("forward_attempt_level.csv","forward_arm_summary.csv","forward_strata.csv",
                             "forward_charge_sensitivity.csv")}}
    (a.output/"forward_evaluation.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
