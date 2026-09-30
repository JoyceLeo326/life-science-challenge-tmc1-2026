"""Outcome-blind checks and file hashes for corrected prospective freeze."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from portable_config import required_path
from rdkit import Chem

HERE=Path(__file__).resolve().parent
ROOT=required_path("TMC1_PREFLIGHT_ROOT")
def read(p):
    with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
arms=read(ROOT/"arm_assignments.csv")
union=read(ROOT/"physical_union_docking_list.csv")
freeze=json.loads((ROOT/"freeze.json").read_text(encoding="utf-8"))
assert len(arms)==300 and len(union)==278
assert Counter(r["arm"] for r in arms)=={"random":120,"descriptor_ridge":120,"AI_initial":60}
assert len({r["compound_id"] for r in union})==278
assert len({r["parent_key"] for r in union})==278
assert all(r["docking_eligibility"]=="eligible" for r in union)
assert all(r["source_stereo_unspecified"]=="False" and r["rule_stable_stereo_unspecified"]=="False" for r in union)
assert all(r["historical_label_overlap"]=="False" for r in union)
charges=Counter(int(float(r["formal_charge"])) for r in union)
rule_charge_by_id={r["compound_id"]:Chem.GetFormalCharge(Chem.MolFromSmiles(r["rule_smiles"])) for r in union}
by_arm={a:{"n":sum(r["arm"]==a for r in arms),
           "exchangeable_N_approx":sum(r["arm"]==a and r["exchangeable_N_geometry_approximation"]=="True" for r in arms),
           "source_abs_charge_ge2":sum(r["arm"]==a and abs(int(float(r["formal_charge"])))>=2 for r in arms),
           "rule_abs_charge_ge2":sum(r["arm"]==a and abs(rule_charge_by_id[r["compound_id"]])>=2 for r in arms)}
        for a in ("random","descriptor_ridge","AI_initial")}
report={"status":"v2 corrected freeze before any v2 docking outcome",
        "source_library_sha256":freeze["input_sha256"],"arm_rows":len(arms),
        "physical_unique":len(union),"selection_overlap":len(arms)-len(union),
        "by_arm":by_arm,"source_formal_charge_counts":dict(charges),
        "rule_formal_charge_counts":dict(Counter(rule_charge_by_id.values())),
        "all_source_and_rule_stable_stereo_specified":True,
        "files_sha256":{p.name:sha(p) for p in [ROOT/"all_predictions.csv",ROOT/"arm_assignments.csv",ROOT/"physical_union_docking_list.csv",ROOT/"freeze.json"]},
        "aborted_prior_freeze":"ranking_formal_16698, 175 preparation logs and zero pose outputs; selection changed before any new docking score existed"}
(ROOT/"preflight_v2.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
print(report)
