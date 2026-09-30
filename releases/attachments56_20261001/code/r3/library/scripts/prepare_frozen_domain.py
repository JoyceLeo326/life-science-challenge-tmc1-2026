"""Recover the exact frozen 9098-ID domain without changing any library or rank."""
import csv,hashlib,json,shutil,sys
from pathlib import Path
from datetime import datetime,timezone,timedelta
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
source=Path(sys.argv[1])
dest=ROOT/'frozen_domain';dest.mkdir(parents=True,exist_ok=True)
expected='79b775967bfc87ded9a07aae03549c9233691ec4fc3375a2c2662b64ee4b2841'
assert sha(source/'all_predictions.csv')==expected
preflight=json.loads((source/'preflight_v2.json').read_text(encoding='utf-8-sig'))
assert preflight['files_sha256']['all_predictions.csv']==expected
rows=list(csv.DictReader((source/'all_predictions.csv').open(encoding='utf-8-sig',newline='')))
eligible=[x for x in rows if x['docking_eligibility']=='eligible']
assert len(rows)==16698 and len(eligible)==9098
assert len({x['compound_id'] for x in eligible})==9098
assert len({x['parent_key'] for x in eligible})==9098
for name in ('all_predictions.csv','preflight_v2.json','freeze.json'):
    shutil.copyfile(source/name,dest/name)
fields=['compound_id','source_id','parent_key','inchi_key','standardized_isomeric_smiles','rule_smiles','rule_inchikey','mw','clogp','tpsa','hbd','hba','rotatable_bonds','formal_charge','docking_eligibility']
with (dest/'eligible_9098.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fields);w.writeheader();w.writerows({k:x[k] for k in fields} for x in eligible)
ids=''.join(x['compound_id']+'\n' for x in eligible)
(dest/'eligible_9098_ids.txt').write_text(ids,encoding='utf-8')
out={'frozen_at_cst':datetime.now(timezone(timedelta(hours=8))).isoformat(),'source':'Frozen R2 ranking_formal_v2_16698 all_predictions.csv, bound by original preflight_v2.json','all_predictions_sha256':expected,'selection_predicate':"docking_eligibility == 'eligible'",'forward_rows':16698,'eligible_rows':9098,'unique_compound_ids':9098,'unique_parent_keys':9098,'excluded_by_reason':{v:sum(x['docking_eligibility']==v for x in rows) for v in sorted({x['docking_eligibility'] for x in rows}) if v!='eligible'},'eligible_csv_sha256':sha(dest/'eligible_9098.csv'),'ordered_ids_sha256':sha(dest/'eligible_9098_ids.txt'),'primary_representation':'rule_smiles: exact frozen pH7.4 docking-rule representation','predeclared_representation_sensitivity':'standardized_isomeric_smiles as standardized source representation; report both separately and union as a diagnostic alert, never as exclusion','changes_to_original_library_or_rank':False}
(dest/'DOMAIN_FREEZE.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print('Frozen domain:',len(eligible),out['eligible_csv_sha256'])
