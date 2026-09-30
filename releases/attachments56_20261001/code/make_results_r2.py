"""Recompute published endpoints from immutable attempt and independent-QC tables.

Python standard library only. This re-evaluates existing calculations; it does
not rerun docking, training, or wet experiments. Paths are package-relative.
"""
import argparse, csv, hashlib, json, time
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parent
THRESHOLD=-8.07960033416748
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def flag(x):
    if x not in ('True','False'):raise ValueError('Nonexplicit boolean')
    return x=='True'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=ROOT/'results.csv');ap.add_argument('--progress-delay',type=float,default=0)
    a=ap.parse_args();start=time.perf_counter()
    def status(s):print(s,flush=True);time.sleep(a.progress_delay)
    index=json.loads((ROOT/'SOURCE_INDEX.json').read_text(encoding='utf-8'))
    for r in index:
        if hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()!=r['sha256']:raise ValueError('Input hash changed: '+r['path'])
    status('Input integrity PASS: '+str(len(index))+' bound files')
    data=rows(ROOT/'data/forward_attempt_level.csv')
    assert len(data)==360 and len({r['compound_id'] for r in data})==319
    qc={r['compound_id']:r for n in ['formal_v2_initial_independent_qc.csv','formal_v2_active_independent_qc.csv'] for r in rows(ROOT/'05_independent_audit'/n)}
    for r in data:
        q=qc[r['compound_id']];assert flag(q['independent_qc_pass']) and q['run_id']==r['run_id']
        assert abs(float(q['score_kcal_mol'])-float(r['score_kcal_mol']))<1e-6
        assert flag(r['favorable_score_predeclared'])==(float(r['score_kcal_mol'])<=THRESHOLD)
        assert flag(r['middle_contact_4A'])==flag(q['middle_contact_4A'])
        assert flag(r['favorable_and_middle_contact'])==(flag(r['favorable_score_predeclared']) and flag(r['middle_contact_4A']))
    assert Counter('AI' if r['arm'].startswith('AI_') else r['arm'] for r in data)=={'random':120,'descriptor_ridge':120,'AI':120}
    status('R2 PASS: 360 attempts / 319 unique physical jobs / 120 per arm')
    results=[]
    for scope in ['R2_all','R3_MW350_matched']:
        for arm in ['random','descriptor_ridge','AI']:
            subset=[r for r in data if ('AI' if r['arm'].startswith('AI_') else r['arm'])==arm and (scope=='R2_all' or float(r['rule_mw'])>=350)]
            results.append({'scope':scope,'arm':arm,'attempted':len(subset),'qc_pass':sum(flag(r['independent_qc_pass']) for r in subset),'favorable_score':sum(flag(r['favorable_score_predeclared']) for r in subset),'middle_contact_4A':sum(flag(r['middle_contact_4A']) for r in subset),'favorable_and_middle_contact':sum(flag(r['favorable_and_middle_contact']) for r in subset)})
    matched={r['arm']:r for r in results if r['scope']=='R3_MW350_matched'}
    assert [(matched[k]['attempted'],matched[k]['favorable_score'],matched[k]['middle_contact_4A']) for k in ['random','descriptor_ridge','AI']]==[(81,35,70),(120,86,119),(120,115,107)]
    status('MW >= 350 PASS: Random 35/81; Ridge 86/120; AI 115/120')
    conf=rows(ROOT/'05_independent_audit/confirmation_four_condition_independent_table.csv')
    assert len(conf)==24 and len({r['compound_id'] for r in conf})==6 and len({r['batch'] for r in conf})==4
    assert sum(flag(r['qc_pass']) for r in conf)==24 and sum(flag(r['middle_contact_4A']) for r in conf)==22
    results.append({'scope':'R2_confirmation_four_conditions','arm':'six_candidates','attempted':24,'qc_pass':24,'favorable_score':sum(float(r['score_kcal_mol'])<=THRESHOLD for r in conf),'middle_contact_4A':22,'favorable_and_middle_contact':22})
    status('Confirmation PASS: 6 candidates x 4 conditions; QC 24/24; Middle 22/24')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    summary_path=a.output.parent/'summary_metrics.csv'
    with summary_path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,results[0].keys());w.writeheader();w.writerows(results)
    panel=rows(ROOT/'data/confirmation_panel.csv');assert len(panel)==6
    assert {r['compound_id'] for r in panel}=={r['compound_id'] for r in conf}
    candidates=[]
    for p in panel:
        cid=p['compound_id'];baseline=[r for r in data if r['compound_id']==cid];assert baseline
        assert all(abs(float(r['score_kcal_mol'])-float(p['initial_score_kcal_mol']))<1e-6 for r in baseline)
        c=[r for r in conf if r['compound_id']==cid];assert len(c)==4
        structure='structures/'+cid+'.sdf';assert (ROOT/structure).is_file()
        candidates.append({'candidate_id':cid,'track':'AI-assisted small-molecule drug design','standardized_isomeric_smiles':p['standardized_isomeric_smiles'],'docking_rule_smiles':p['rule_smiles'],'source_id':p['source_id'],'selection_rank':p['confirmation_selection_rank'],'initial_vina_kcal_mol':p['initial_score_kcal_mol'],'ExtraTrees_predicted_vina_kcal_mol':p['ai_pred_kcal_mol'],'tree_SD_kcal_mol':p['ai_tree_sd_kcal_mol'],'confirmation_QC_pass':sum(flag(r['qc_pass']) for r in c),'confirmation_conditions':4,'confirmation_middle_contact_pass':sum(flag(r['middle_contact_4A']) for r in c),'confirmation_scores_kcal_mol':';'.join(r['batch']+':'+r['score_kcal_mol'] for r in c),'model_version':'ExtraTrees160_seed20260929_Morgan1024_8descriptors','run_version':'R2_M_PUB_ALL_Vina1.2.7_e8_seed20260927;confirmation_four_conditions','structure_file':structure,'remarks':'Frozen three-arm favorable/QC/Middle panel; one per Murcko scaffold. Static computational evidence; TMC1 functional direction and USH1B compensation remain untested.'})
    with a.output.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,candidates[0].keys());w.writeheader();w.writerows(candidates)
    generated=rows(ROOT/'data/final_generated40_actual_outcomes.csv');assert len(generated)==40 and sum(flag(r['joint_qc_score_middle_pass']) for r in generated)==0
    status('Candidate list PASS: 6 verified identities, SMILES, versions and ligand structures; original generation joint 0/40 retained')
    report={'status':'passed','operation':'re-evaluation of existing scores and independent-QC geometry; no new docking','threshold':THRESHOLD,'input_files':index,'results_sha256':hashlib.sha256(a.output.read_bytes()).hexdigest(),'wall_seconds':time.perf_counter()-start,'checks':{'attempts_360':True,'unique_physical_319':True,'per_arm_120':True,'QC_score_geometry_join':True,'MW350_matched':True,'confirmation_24':True}}
    (a.output.parent/'RUN_REPORT.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    status('SUCCESS: 6-candidate results.csv and separate summary_metrics.csv written; SHA256 '+report['results_sha256'])
if __name__=='__main__':main()
