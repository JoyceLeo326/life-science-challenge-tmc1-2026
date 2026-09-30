"""Read-only acceptance of the adjacent A1-A3 evidence directories; needs RDKit."""
import csv,json,hashlib,math
from pathlib import Path
from rdkit import Chem
from rdkit.Contrib.SA_Score import sascorer
B=Path(__file__).resolve().parent
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
checks=[]
def check(name,v):
    checks.append({'check':name,'pass':bool(v)})
    if not v:raise AssertionError(name)
for r in json.loads((B/'FILE_MANIFEST.json').read_text(encoding='utf-8'))['files']:
    p=B/r['path'];check('SHA256:'+r['path'],p.is_file() and p.stat().st_size==r['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'])
a1=B/'A1_reference_sources';a2=B/'A2_chembl_neighbors';a3=B/'A3_ADMET_SA'
panel=read(a3/'ADMET_SA_panel11.csv');sw={r['Molecule']:r for r in read(a3/'SwissADME_combined11.csv')};sa={r['id']:r for r in read(a3/'SA_physchem_actual11.csv')}
check('11 unique SwissADME IDs',len(panel)==len(sw)==11)
for r in panel:
    check('SwissADME identity '+r['id'],Chem.MolToInchiKey(Chem.MolFromSmiles(sw[r['id']]['Canonical SMILES']))==r['inchikey'])
    check('SA actual recompute '+r['id'],abs(sascorer.calculateScore(Chem.MolFromSmiles(r['smiles']))-float(sa[r['id']]['SA_score']))<1e-9)
pairs=read(a2/'candidate_neighbor_annotations.csv');check('1458 pairs',len(pairs)==1458);check('847 unique full keys',len({r['neighbor_inchikey'] for r in pairs})==847)
for r in read(a2/'candidate_similarity_summary.csv'):
    z=[p for p in pairs if p['candidate_chembl_id']==r['candidate_chembl_id']]
    check('40pct count '+r['candidate_chembl_id'],len(z)==int(r['count_ge_0p4_excluding_self']))
    check('60pct count '+r['candidate_chembl_id'],sum(float(p['official_tanimoto'])>=.6 for p in z)==int(r['count_ge_0p6_excluding_self']))
check('622 activities',len(read(a2/'pharmacology_actual_activity_records.csv'))==622)
ai=read(a3/'ADMET_AI_v2_11_full_predictions.csv');meta={r['id']:r for r in read(a3/'ADMET_AI_v2_endpoint_metadata.csv')};exc=0;n=0
for r in ai:
    for k,v in r.items():
        if k not in meta:continue
        x=float(v);m=meta[k];check(r['id']+' finite '+k,math.isfinite(x));n+=1
        exc+=x<float(m['minimum']) or x>float(m['maximum'])
check('451 predicted endpoints',n==451);check('10 raw range exceptions',exc==10 and len(read(a3/'ADMET_AI_v2_out_of_range.csv'))==10)
check('unverified68 not converted into members',read(a1/'reference_source_claim_audit.csv')[0]['status']=='SOURCE_NOT_ESTABLISHED')
print(json.dumps({'passed':len(checks),'failed':0,'checks':checks},ensure_ascii=False))
