"""Offline reconstruction from saved API responses and actual prediction CSVs."""
import argparse,csv,json,hashlib
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Contrib.SA_Score import sascorer
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def jsons(folder,pattern,key):
    rows=[]
    for p in sorted(folder.glob(pattern)):
        if p.name.endswith('.receipt.json'):continue
        rows+=json.loads(p.read_text(encoding='utf-8'))[key]
    return rows
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args();out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    b=Path(__file__).resolve().parent;a1=b/'A1_reference_sources';a2=b/'A2_chembl_neighbors';a3=b/'A3_ADMET_SA';f=a2/'pharmacology_raw'
    acts=jsons(f,'activity_batch*_page*.json','activities');docs={r['document_chembl_id']:r for r in jsons(f,'document_batch*_page*.json','documents')};assays={r['assay_chembl_id']:r for r in jsons(f,'assay_batch*_page*.json','assays')}
    activity=[]
    for a in acts:
        d=docs.get(a.get('document_chembl_id'),{});s=assays.get(a.get('assay_chembl_id'),{});keys=['activity_id','molecule_chembl_id','molecule_pref_name','target_chembl_id','target_pref_name','target_organism','assay_chembl_id','assay_type','assay_description','standard_type','standard_relation','standard_value','standard_units','pchembl_value','standard_flag','data_validity_comment','potential_duplicate','document_chembl_id']
        r={k:a.get(k) for k in keys};r.update(assay_confidence_score=s.get('confidence_score'),assay_bao_format=s.get('bao_format'),document_doi=d.get('doi'),document_title=d.get('title'),document_year=d.get('year'),document_journal=d.get('journal'),source_api=f'https://www.ebi.ac.uk/chembl/api/data/activity/{a["activity_id"]}.json');activity.append(r)
    write(out/'pharmacology_actual_activity_records.csv',activity)
    pairs=[]
    for p in read(a2/'formal_six_input.csv'):
        cid=p['source_id'];rr=jsons(a2/'chembl_raw'/cid,'similarity40_page*.json','molecules')
        for r in sorted(rr,key=lambda r:(-float(r['similarity']),r['molecule_chembl_id'])):
            st=r['molecule_structures'];pairs.append(dict(candidate_id=p['id'],candidate_chembl_id=cid,neighbor_chembl_id=r['molecule_chembl_id'],official_tanimoto=float(r['similarity'])/100,threshold_ge_0p4=float(r['similarity'])>=40,threshold_ge_0p6=float(r['similarity'])>=60,neighbor_inchikey=st['standard_inchi_key'],neighbor_smiles=st['canonical_smiles']))
    write(out/'two_threshold_neighbors_rebuilt.csv',pairs)
    sw={r['Molecule']:r for name in ['swissadme_actual.csv','swissadme_remaining6_actual.csv'] for r in read(a3/name)};ai={r['id']:r for r in read(a3/'ADMET_AI_v2_11_full_predictions.csv')};rows=[]
    panel=read(a3/'ADMET_SA_panel11.csv')
    for r in panel:
        s=sw[r['id']];assert Chem.MolToInchiKey(Chem.MolFromSmiles(s['Canonical SMILES']))==r['inchikey']
        q={x:r[x] for x in ('id','name','group','pair_role','inchikey')};q['RDKit_SA_score']=sascorer.calculateScore(Chem.MolFromSmiles(r['smiles']))
        for x in ('MW','Consensus Log P','TPSA','GI absorption','BBB permeant','Pgp substrate','CYP1A2 inhibitor','CYP2C19 inhibitor','CYP2C9 inhibitor','CYP2D6 inhibitor','CYP3A4 inhibitor','ESOL Log S','ESOL Class','Lipinski #violations','PAINS #alerts','Brenk #alerts','Synthetic Accessibility'):q['SwissADME_'+x]=s[x]
        for x in ('AMES','hERG','DILI','ClinTox','CYP3A4_Veith','Solubility_AqSolDB','PPBR_AZ','VDss_Lombardo'):q['ADMET_AI_'+x]=ai[r['id']][x]
        rows.append(q)
    write(out/'multiparameter_druglikeness_rebuilt.csv',rows)
    def par(s):return Chem.MolToInchiKey(rdMolStandardize.Uncharger().uncharge(rdMolStandardize.FragmentParent(Chem.MolFromSmiles(s))))
    keys={par(r['representative_isomeric_smiles']) for r in read(a1/'author_training_unique12.csv')};overlap=[dict(candidate_id=r['id'],candidate_parent_inchikey=par(r['smiles']),author_training12_parent_overlap=par(r['smiles']) in keys) for r in panel if r['group']=='formal_six'];write(out/'formal6_training12_overlap_rebuilt.csv',overlap)
    # Match preserved complete source table and every shared field of the reduced reconstruction.
    comparisons=[]
    for new,old,idfields in [('pharmacology_actual_activity_records.csv',a2/'pharmacology_actual_activity_records.csv',['activity_id']),('two_threshold_neighbors_rebuilt.csv',a2/'candidate_neighbor_annotations.csv',['candidate_id','neighbor_chembl_id']),('multiparameter_druglikeness_rebuilt.csv',a3/'multiparameter_druglikeness11.csv',['id']),('formal6_training12_overlap_rebuilt.csv',a1/'formal6_author_training12_overlap.csv',['candidate_id'])]:
        aa=read(out/new);bb=read(old);bb={tuple(r[k] for k in idfields):r for r in bb};differences=[]
        for r in aa:
            oldr=bb[tuple(r[k] for k in idfields)]
            for k,v in r.items():
                if v!=oldr[k]:
                    try:ok=abs(float(v)-float(oldr[k]))<1e-9
                    except ValueError:ok=False
                    if not ok:differences.append([r.get('id',r.get('activity_id',r.get('candidate_id'))),k,v,oldr[k]])
        assert not differences and len(aa)==len(bb),(new,differences[:3]);comparisons.append(dict(table=new,rows=len(aa),differing_shared_fields=0))
    receipt=dict(status='PASS',mode='offline_saved_primary_responses_and_predictions',comparisons=comparisons,not_new_docking_or_ADMET_model_training=True)
    (out/'OFFLINE_REBUILD_RECEIPT.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps(receipt))
if __name__=='__main__':main()
