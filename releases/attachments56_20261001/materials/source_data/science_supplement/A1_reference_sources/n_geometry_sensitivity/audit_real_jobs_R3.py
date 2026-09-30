"""Independent read-only chemistry, raw-score, hash, geometry and region audit.

Usage: python audit_real_jobs.py BATCH OUTPUT_STEM
One CPU, no docking. Results are a snapshot of the ledger at read time.
"""
import csv,hashlib,json,math,re,sys
from collections import Counter
from pathlib import Path
import numpy as np
from rdkit import Chem

BASE=Path(r"RAW_RUN_ROOT")
if len(sys.argv)>3:BASE=Path(sys.argv[3])
OUT=Path(__file__).resolve().parent
MIDDLE={408,411,412,443,444,447,448,528,531,532,536,579,580,582,601}
MIDDLE_H={414,417,418,449,450,453,454,531,534,535,539,582,583,585,604}
BATCH=sys.argv[1];STEM=sys.argv[2]
ROOT=BASE/BATCH
EVIDENCE=str(OUT/(STEM+".csv"))
def read(p):
    with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def sdf(p):
    with Path(p).open("rb") as f:return next(Chem.ForwardSDMolSupplier(f,removeHs=False))
def stable(m,indices):
    if m is None:return None
    m=Chem.Mol(Chem.RemoveHs(m))
    if any(i>=m.GetNumAtoms() for i in indices):return None
    for i in indices:m.GetAtomWithIdx(i).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    Chem.AssignStereochemistry(m,cleanIt=True,force=True)
    return Chem.MolToSmiles(m,isomericSmiles=True)
def atom_lines(p,first=False):
    out=[]
    for l in Path(p).read_text(encoding="utf-8",errors="replace").splitlines():
        if first and l.startswith("ENDMDL"):break
        if l.startswith(("ATOM  ","HETATM")):
            t=l.split()[-1]
            if t in {"H","HD","HS"} or re.fullmatch(r"G\d+",t):continue
            try:out.append((l[21],int(l[22:26]),l[17:20].strip(),np.array([float(l[i:i+8]) for i in (30,38,46)])))
            except ValueError:pass
    return out
def log_score(p):
    for l in Path(p).read_text(encoding="latin-1").splitlines():
        s=l.split()
        if len(s)>=4 and s[0]=="1":
            try:
                value=float(s[1])
                if math.isfinite(value):return value
            except ValueError:pass
    return None
def pose_header_score(p):
    for line in Path(p).read_text(encoding="utf-8",errors="replace").splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            try:
                value=float(line.split()[3])
                if math.isfinite(value):return value
            except (ValueError,IndexError):pass
        if line.startswith("ATOM  "):break
    return None
def pdbqt_charge_sum(p):
    charges=[]
    for line in Path(p).read_text(encoding="utf-8",errors="replace").splitlines():
        if line.startswith(("ATOM  ","HETATM")):
            try:charges.append(float(line.split()[-2]))
            except (ValueError,IndexError):return None
    return sum(charges) if charges else None

prep={r["compound_id"]:r for r in read(ROOT/"preparation_ledger.csv")}
jobs=read(ROOT/"job_ledger.csv")
cache={};out=[]
for j in jobs:
    cid=j["compound_id"];p=prep.get(cid,{});bad=[]
    rule=Chem.MolFromSmiles(j.get("rule_smiles",j.get("standardized_isomeric_smiles","")))
    ni=json.loads(j.get("exchangeable_protonated_N_indices") or "[]")
    if p.get("preparation_status")!="prepared":bad.append("preparation_not_prepared")
    if rule is None:bad.append("invalid_rule_smiles")
    if not Path(j.get("pose_file","")).is_file():bad.append("missing_pose")
    if not Path(j.get("log_file","")).is_file():bad.append("missing_vina_log")
    if not Path(j.get("exported_sdf_file","")).is_file():bad.append("missing_export_sdf")
    if bad:
        out.append({"compound_id":cid,"run_id":j.get("run_id",""),"independent_qc_pass":False,
                    "failure_reason":";".join(bad),"score_kcal_mol":"","pose_sha256":"","status":j.get("status",""),
                    "independent_qc_evidence":EVIDENCE});continue
    inputmol=sdf(j["sdf_file"]);outmol=sdf(j["exported_sdf_file"])
    if stable(inputmol,ni)!=stable(rule,ni):bad.append("input_sdf_stable_graph_mismatch")
    if stable(outmol,ni)!=stable(rule,ni):bad.append("pose_sdf_stable_graph_mismatch")
    ligqt=Path(j["ligand_pdbqt"]).read_text(encoding="utf-8",errors="replace")
    charge_sum=pdbqt_charge_sum(j["ligand_pdbqt"])
    charge_delta=abs(charge_sum-Chem.GetFormalCharge(rule)) if charge_sum is not None else float("inf")
    # Vina scoring ignores PDBQT user partial charges (official FAQ). Record
    # charge-conservation anomalies for chemical-state review, not score QC.
    charge_warning=charge_delta>0.05
    remark=next((l.removeprefix("REMARK SMILES ").strip() for l in ligqt.splitlines() if l.startswith("REMARK SMILES ")),"")
    remarkmol=Chem.MolFromSmiles(remark) if remark else None
    if stable(remarkmol,ni)!=stable(rule,ni):bad.append("ligand_pdbqt_remark_graph_mismatch")
    if sha(j["ligand_pdbqt"])!=j["ligand_sha256"]:bad.append("ligand_pdbqt_hash_mismatch")
    if sha(j["sdf_file"])!=j["sdf_sha256"]:bad.append("input_sdf_hash_mismatch")
    if sha(j["receptor_file"])!=j["receptor_sha256"]:bad.append("receptor_hash_mismatch")
    if sha(j["box_file"])!=j["box_sha256"]:bad.append("box_hash_mismatch")
    actual_pose_sha=sha(j["pose_file"])
    if actual_pose_sha!=j["pose_sha256"]:bad.append("pose_hash_mismatch")
    if sha(j["exported_sdf_file"])!=j["exported_sdf_sha256"]:bad.append("pose_sdf_hash_mismatch")
    if j["receptor_file"] not in cache:
        rec=atom_lines(j["receptor_file"])
        site=MIDDLE_H if j.get("receptor_id")=="H_AF" else MIDDLE
        cache[j["receptor_file"]]=(np.array([a[-1] for a in rec]),np.array([a[-1] for a in rec if a[0]=="A" and a[1] in site]))
    rxyz,middle=cache[j["receptor_file"]]
    pose=atom_lines(j["pose_file"],first=True);xyz=np.array([a[-1] for a in pose])
    if len(xyz)!=rule.GetNumHeavyAtoms():bad.append("pose_heavy_atom_count_mismatch")
    mind=float("inf");severe=0
    for i in range(0,len(rxyz),1000):
        d=np.linalg.norm(xyz[:,None,:]-rxyz[None,i:i+1000,:],axis=2)
        mind=min(mind,float(d.min()));severe+=int(np.count_nonzero(d<1.5))
    mdmiddle=float(np.linalg.norm(xyz[:,None,:]-middle[None,:,:],axis=2).min())
    center=np.array(json.loads(j["box_center_A"]));half=np.array(json.loads(j["box_size_A"]))/2
    margin=float(np.min(half-np.abs(xyz-center)))
    if mind<1.5:bad.append("receptor_heavy_clash")
    if margin<-.1:bad.append("outside_docking_box")
    score=pose_header_score(j["pose_file"])
    table_score=log_score(j["log_file"])
    if score is None or not j.get("score_kcal_mol") or abs(score-float(j["score_kcal_mol"]))>1e-6:
        bad.append("pose_header_score_mismatch")
    if table_score is None or score is None or abs(table_score-score)>0.0051:
        bad.append("vina_log_table_score_mismatch")
    if j["status"]!="completed":bad.append("job_not_completed")
    out.append({"compound_id":cid,"run_id":j["run_id"],"independent_qc_pass":not bad,
                "failure_reason":";".join(bad),"score_kcal_mol":score if score is not None else "",
                "pose_sha256":actual_pose_sha,"status":j["status"],
                "independent_qc_evidence":EVIDENCE,
                "rule_charge":Chem.GetFormalCharge(rule),"exchangeable_N_indices":json.dumps(ni),
                "pdbqt_partial_charge_sum":round(charge_sum,3) if charge_sum is not None else "",
                "pdbqt_charge_abs_difference":round(charge_delta,3),
                "charge_parameterization_warning":charge_warning,
                "vina_log_table_score_kcal_mol":table_score if table_score is not None else "",
                "pose_heavy_atoms":len(xyz),"min_receptor_A":round(mind,3),"severe_pairs_under1p5_A":severe,
                "min_box_margin_A":round(margin,3),"closest_middle_A":round(mdmiddle,3),
                "middle_contact_4A":mdmiddle<4,"log_sha256":sha(j["log_file"])})
with (OUT/(STEM+".csv")).open("w",encoding="utf-8-sig",newline="") as f:
    fields=list(dict.fromkeys(k for r in out for k in r));w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
summary={"batch":BATCH,"jobs_in_ledger":len(jobs),"independent_qc_pass":sum(r["independent_qc_pass"] for r in out),
         "independent_qc_failed":sum(not r["independent_qc_pass"] for r in out),
         "reasons":dict(Counter(reason for r in out for reason in r["failure_reason"].split(";") if reason)),
         "middle_contact_4A_among_qc_pass":sum(r["independent_qc_pass"] and r.get("middle_contact_4A",False) for r in out),
         "charge_parameterization_warnings":sum(r.get("charge_parameterization_warning",False) for r in out),
         "score_le_frozen_threshold_among_qc_pass":sum(r["independent_qc_pass"] and r.get("score_kcal_mol","")!="" and float(r["score_kcal_mol"])<=-8.07960033416748 for r in out),
         "source_ledger_file_sha256":sha(ROOT/"job_ledger.csv")}
(OUT/(STEM+"_summary.json")).write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps(summary,ensure_ascii=False,indent=2))
