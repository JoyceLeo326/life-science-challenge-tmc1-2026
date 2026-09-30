"""Independent full-record check of 02's exact-stereo SDF and exclusions."""
import csv,hashlib,json
from pathlib import Path
import os
from rdkit import Chem

ROOT=Path(os.environ["TMC1_R2_SOURCE_ROOT"])
CSV=ROOT/"02_公共分子库/processed/filtered_for_screening.csv"
SNAP=ROOT/"02_公共分子库/snapshots/v2_20000_new"
SDF=SNAP/"filtered_for_screening_exact_stereo_fullprops_2d.sdf"
EXC=SNAP/"sdf_exact_stereo_exclusions.csv"
MAN=SNAP/"sdf_exact_stereo_manifest.json"
OUT=ROOT/"05_独立方法复核"
def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
source=rows(CSV);exclusions=rows(EXC);manifest=json.loads(MAN.read_text(encoding="utf-8"))
excluded={r["compound_id"] for r in exclusions}
remaining=[r for r in source if r["compound_id"] not in excluded]
prediction={r["compound_id"]:r for r in rows(ROOT/"03_算法与计算/ranking_formal_v2_16698/all_predictions.csv")}
bad=[];count=0;property_names=[]
with SDF.open("rb") as f:
    for count,mol in enumerate(Chem.ForwardSDMolSupplier(f,removeHs=True),start=1):
        if count>len(remaining):bad.append([count,"extra_record"]);continue
        r=remaining[count-1]
        if mol is None:bad.append([count,r["compound_id"],"parse_failed"]);continue
        if count==1:property_names=list(mol.GetPropNames())
        if Chem.MolToInchiKey(mol)!=r["inchi_key"]:
            bad.append([count,r["compound_id"],"inchi_key_mismatch"])
        source_mol=Chem.MolFromSmiles(r["standardized_isomeric_smiles"])
        if Chem.MolToSmiles(mol,isomericSmiles=False)!=Chem.MolToSmiles(source_mol,isomericSmiles=False):
            bad.append([count,r["compound_id"],"connectivity_mismatch"])
        for field in ("compound_id","source_id","source_ids","raw_smiles",
                      "standardized_isomeric_smiles","inchi_key","parent_key",
                      "source_url","processing_status"):
            if mol.HasProp(field) and mol.GetProp(field)!=r[field]:
                bad.append([count,r["compound_id"],f"property_{field}_mismatch"])
previous=json.loads((OUT/"fullprops_sdf_graphs_independent_audit.json").read_text(encoding="utf-8"))
previous_bad={r["compound_id"] for r in previous["all_key_mismatches"]}
result={"source_csv_rows":len(source),"excluded_rows":len(exclusions),"excluded_unique":len(excluded),
        "repaired_sdf_records":count,"expected_remaining":len(remaining),
        "all_remaining_exact_keys_and_connectivity_match":count==len(remaining) and not bad,
        "record_or_property_errors":len(bad),"first_errors":bad[:10],"first_record_properties":property_names,
        "exclusions_subset_of_prior_key_mismatches":excluded<=previous_bad,
        "all_exclusions_ineligible_for_formal_docking":all(prediction[cid]["docking_eligibility"]=="unspecified_source_stereo" for cid in excluded),
        "prior_key_mismatches_repaired":len(previous_bad-excluded),
        "manifest_hashes_match":sha(CSV)==manifest["source_csv_sha256"] and sha(SDF)==manifest["repaired_sdf_sha256"] and sha(EXC)==manifest["excluded_csv_sha256"],
        "repaired_sdf_sha256":sha(SDF),"exclusions_sha256":sha(EXC)}
(OUT/"repaired_sdf_independent_audit.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps(result,ensure_ascii=False,indent=2))
if not (result["all_remaining_exact_keys_and_connectivity_match"] and result["exclusions_subset_of_prior_key_mismatches"] and result["all_exclusions_ineligible_for_formal_docking"] and result["manifest_hashes_match"]):raise SystemExit(1)
