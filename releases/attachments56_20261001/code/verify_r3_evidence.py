"""Recompute R3 summaries from frozen evidence. Python standard library only."""
import argparse,csv,json,hashlib
from pathlib import Path
from collections import Counter

ROOT=Path(__file__).resolve().parent
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def yes(x):return str(x).lower() in ('true','1')
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=ROOT/'r3_summary.csv');a=ap.parse_args()
    snapshot=json.loads((ROOT/'INTEGRATED_SOURCE_SNAPSHOT_INDEX.json').read_text(encoding='utf-8'))
    bound=[r for r in snapshot['files'] if r['component']!='R1_R2_formal_code' and r.get('included',True)]
    for r in bound:
        assert sha(ROOT/r['integrated_relative_path'])==r.get('integrated_sha256',r['copied_sha256']),r['integrated_relative_path']
    candidates={r['candidate_id'] for r in rows(ROOT/'results.csv')};assert len(candidates)==6
    library=ROOT/'r3/library';diagnostic=rows(library/'diagnostics/chemical_alerts_9098.csv');domain=rows(library/'frozen_domain/eligible_9098.csv')
    assert len(diagnostic)==len(domain)==9098
    assert len({r['compound_id'] for r in diagnostic})==len({r['parent_key'] for r in diagnostic})==9098
    assert {r['compound_id'] for r in diagnostic}=={r['compound_id'] for r in domain}
    summary=[]
    claimed=rows(library/'diagnostics/alert_summary.csv')
    claim={(r['representation'],r['family']):r for r in claimed}
    families=('pains','pains_A','pains_B','pains_C','brenk','reactive')
    for rep in ('primary','sensitivity','union'):
        for family in families:
            counts=[]
            for r in diagnostic:
                rules=json.loads(r[f'{rep}_{family}_rules_json']);assert len(rules)==int(r[f'{rep}_{family}_rule_count'])
                if rep=='union':assert set(rules)==set(json.loads(r[f'primary_{family}_rules_json']))|set(json.loads(r[f'sensitivity_{family}_rules_json']))
                counts.append(bool(rules))
            hits=sum(counts);c=claim[(rep,family)]
            assert hits==int(c['hit_molecules']) and int(c['denominator_all_attempted'])==9098
            failed=sum(r['overall_status']!='OK' if rep=='union' else r[rep+'_status']!='OK' for r in diagnostic)
            assert failed==int(c['failed_or_incomplete_molecules'])==0
            summary.append({'component':'chemical_alerts','metric':rep+'_'+family+'_hit_molecules','value':hits,'denominator':9098,'scope':'structural alert, not biological activity or measured toxicity'})
    six=[r for r in diagnostic if r['compound_id'] in candidates];assert len(six)==6
    assert {r['compound_id'] for r in six}=={r['compound_id'] for r in diagnostic if yes(r['is_six_candidate'])}
    assert all(int(r[f'{rep}_{family}_rule_count'])==0 for r in six for rep in ('primary','sensitivity','union') for family in families)
    matrix=rows(ROOT/'r3/structure/six_candidates_plus_two_references_three_species24.csv')
    jobs=rows(ROOT/'r3/structure/all26_jobs_anonymous.csv')
    assert len(matrix)==24 and len({(r['compound_id'],r['structure_id']) for r in matrix})==24
    assert len({r['compound_id'] for r in matrix})==8 and len({r['structure_id'] for r in matrix})==3
    assert sum(yes(r['independent_qc_pass']) for r in matrix)==23
    assert len(jobs)==len({r['run_id'] for r in jobs})==26 and sum(yes(r['independent_qc_pass']) for r in jobs)==25
    worm=[r for r in matrix if r['structure_id']=='W_7USX_LOCAL_COMPLETE'];assert len(worm)==8
    assert all(yes(r['omitted_context_clash_under1p5A']) for r in worm)
    assert {r['compound_id'] for r in matrix if r['row_role']=='candidate'}==candidates
    for metric,value,den in [('static_matrix_protein_QC',23,24),('all_raw_jobs_protein_QC',25,26),('worm_matrix_omitted_HET_clash',8,8)]:
        summary.append({'component':'structure','metric':metric,'value':value,'denominator':den,'scope':'protein-only technical QC; worm omitted-HET conflict retained; not cross-species affinity'})
    networks=rows(ROOT/'r3/anm/network_numerical_QC.csv');frames=rows(ROOT/'r3/anm/conformation_diagnostics.csv')
    assert len(networks)==6 and all(yes(r['ANM_numerical_QC_pass']) and int(r['zero_mode_count'])==6 for r in networks)
    assert len(frames)==24 and not any(yes(r['eligible_for_full_atom_docking']) for r in frames)
    for amplitude,count in [(0.25,8),(0.5,6),(1.0,0)]:
        f=[r for r in frames if float(r['declared_CA_RMS_displacement_A'])==amplitude];assert len(f)==8 and sum(yes(r['coarse_reference_diagnostic_pass']) for r in f)==count
        summary.append({'component':'ANM','metric':str(amplitude)+'A_coarse_geometry_pass','value':count,'denominator':8,'scope':'C-alpha deterministic traversal only; no MD and no full-atom docking receptor'})
    top=rows(ROOT/'r3/reference_generation/R3_top10_replay_table.csv');assert len(top)==10
    assert sum(r['preparation_status']=='blocked' for r in top)==1
    assert sum(yes(r['independent_qc_pass']) for r in top)==9 and sum(yes(r['joint_pass']) for r in top)==0
    summary.append({'component':'generation_supplement','metric':'joint_pass_selected','value':0,'denominator':10,'scope':'frozen cached proxy top10;9 docked,1 identity blocked; original40 unchanged'})
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,list(summary[0]));w.writeheader();w.writerows(summary)
    report={'status':'PASS','bound_R3_source_files_verified':len(bound),'diagnostic_rows':9098,'all_failures_retained':True,'six_candidates_alert_free':True,'static_matrix':24,'raw_jobs':26,'ANM_networks':6,'ANM_traversals':24,'full_atom_dynamic_receptors':0,'generation_top10':10,'results_csv_sha256':sha(ROOT/'results.csv'),'r3_summary_sha256':sha(a.output),'scope':'existing frozen calculation readback; no new docking/training/MD'}
    a.output.with_name('R3_READBACK_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False))
if __name__=='__main__':main()
