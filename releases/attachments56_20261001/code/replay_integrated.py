"""Run the complete existing R2/R3 replay in an independent output directory.

Frozen data and models remain unchanged. External author receptor files are
required by SHA; neither an installation nor a full docking rerun is performed.
"""
import argparse,csv,hashlib,json,os,shutil,subprocess,sys,time,importlib.util
from pathlib import Path
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def table(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',required=True,type=Path);p.add_argument('--work',required=True,type=Path)
    p.add_argument('--skip-full-ranking',action='store_true');p.add_argument('--skip-small-training',action='store_true');p.add_argument('--resume',action='store_true',help='Continue verified earlier steps from this same independent replay work directory');a=p.parse_args()
    required=['rdkit','numpy','scipy','sklearn','joblib','threadpoolctl','dimorphite_dl','meeko','matplotlib']
    if not a.skip_small_training:required+=['torch','selfies']
    missing=[name for name in required if importlib.util.find_spec(name) is None]
    if missing:raise RuntimeError('Missing replay dependencies: '+', '.join(missing)+'. Use a separate environment with requirements.txt installed.')
    work=a.work.resolve()
    if work==ROOT or work.is_relative_to(ROOT) or ROOT.is_relative_to(work):raise ValueError('Replay work must be separate from the package')
    if work.exists() and any(work.iterdir()) and not a.resume:raise FileExistsError('New or empty replay work required; explicit --resume for the same run')
    inputs=json.loads(a.inputs.read_text(encoding='utf-8'))
    for k in ('mouse_pdbqt','human_pdbqt','mouse_aligned_pdb','human_aligned_pdb'):
        value=Path(inputs[k]).expanduser()
        if not value.is_absolute():value=a.inputs.resolve().parent/value
        if not value.is_file():raise FileNotFoundError('Required external input: '+k)
        inputs[k]=str(value.resolve())
    work.mkdir(parents=True,exist_ok=True);(work/'logs').mkdir(exist_ok=True)
    receipts=json.loads((work/'ACTUAL_REPLAY_RECEIPTS.json').read_text(encoding='utf-8')) if a.resume and (work/'ACTUAL_REPLAY_RECEIPTS.json').is_file() else []
    completed={r['id'] for r in receipts if r['exit_code']==0};checks=[]
    env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONUTF8='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    def public_command(cmd):
        pairs=sorted([(str(ROOT),'$PACKAGE'),(str(work),'$WORK')]+[(v,'$EXTERNAL_'+k.upper()) for k,v in inputs.items()],key=lambda x:-len(x[0]))
        out=[]
        for s in map(str,cmd):
            if s==sys.executable:s='python'
            for value,label in pairs:s=s.replace(value,label)
            out.append(s)
        return out
    def run(name,script,args=(),cwd=ROOT,extra=None):
        if name in completed:print(name,'earlier successful step retained; outputs checked below',flush=True);return
        cmd=[sys.executable,'-B',str(script)]+list(map(str,args));started=datetime.now(timezone.utc).isoformat();t=time.perf_counter()
        runenv=env.copy();runenv.update(extra or {})
        log=work/'logs'/f'{name}.log'
        with log.open('w',encoding='utf-8') as f:r=subprocess.run(cmd,cwd=cwd,env=runenv,stdout=f,stderr=subprocess.STDOUT)
        receipt={'id':name,'command':public_command(cmd),'script_sha256':sha(script),'started_utc':started,'elapsed_seconds':round(time.perf_counter()-t,3),'exit_code':r.returncode,'log_sha256':sha(log)}
        receipts.append(receipt);(work/'ACTUAL_REPLAY_RECEIPTS.json').write_text(json.dumps(receipts,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        print(name,'exit',r.returncode,'seconds',receipt['elapsed_seconds'],flush=True)
        if r.returncode:raise RuntimeError('Failed replay step '+name+'; inspect independent output logs')
    def same(name,source,result):
        ok=sha(source)==sha(result);checks.append({'check':name,'pass':ok,'source_sha256':sha(source),'replay_sha256':sha(result)});assert ok,name
    def scientific_table_same(name,source,result,metadata_fields):
        left,right=table(source),table(result);assert len(left)==len(right)
        fields=set(left[0])-set(metadata_fields);assert fields==set(right[0])-set(metadata_fields)
        ok=all(all(l[k]==r[k] for k in fields) for l,r in zip(left,right))
        checks.append({'check':name,'pass':ok,'source_sha256':sha(source),'replay_sha256':sha(result),'comparison':'every scientific cell exact; published evidence-path metadata normalized','excluded_metadata_fields':metadata_fields});assert ok,name
    run('official_candidates',ROOT/'make_results.py',['--output',work/'candidate/results.csv'])
    same('six_candidate_CSV_byte_identical',ROOT/'results.csv',work/'candidate/results.csv')
    same('R2_summary_byte_identical',ROOT/'summary_metrics.csv',work/'candidate/summary_metrics.csv')
    bundle=ROOT/'reproduction_bundle';r1=bundle/'science/R1';r2=bundle/'science/R2/evidence';prepare=bundle/'reproduction/prepare_portable.py';portable=work/'portable'
    prep_args=['--r1',r1,'--r2',r2,'--work',portable,'--receptor',inputs['mouse_pdbqt']]
    run('prepare_check',prepare,prep_args+['--check-only']);run('prepare_actual',prepare,prep_args)
    p03=portable/'portable_method/03_compute';raw=r2/'03_compute/raw_runs';freeze=r2/'03_compute/ranking_formal_v2_16698';active=r2/'03_compute/ranking_active_v2_after_ai60'
    if not a.skip_full_ranking:
        run('rank_all16698',p03/'rank_library.py',['--library',r2/'02_library/formal_snapshot/forward_unlabeled_pool.csv','--output',work/'ranking_replay','--per-arm','120','--ai-initial','60','--exclude-ids',p03/'runtime_pilot_12.csv','--expected-library-sha256','8fd8424797e441a3921ba3805dd8e9bf87d3421e5fe3440a32211753ec964c1e'])
        for n in ('all_predictions.csv','arm_assignments.csv','physical_union_docking_list.csv'):same('R2_ranking_'+n,freeze/n,work/'ranking_replay'/n)
        run('rank_preflight',p03/'freeze_preflight_v2.py',extra={'TMC1_PREFLIGHT_ROOT':str(work/'ranking_replay')})
    run('active_acquire60',p03/'active_acquire.py',['--freeze',freeze,'--ai-attempts',r2/'05_independent_audit/AI_initial_60_independent_outcomes.csv','--output',work/'active_replay'])
    for n in ('active_predictions_all_eligible.csv','AI_active_60.csv'):same('R2_active_'+n,active/n,work/'active_replay'/n)
    run('active_physical',p03/'prepare_active_physical.py',['--freeze',freeze,'--active',active,'--output',work/'active_physical'])
    run('initial_merge',p03/'merge_formal_batches.py',['--output',work/'merged_initial','--require-complete'])
    run('active_resume',p03/'freeze_active_resume.py',['--output',work/'active_resume'])
    run('active_merge',p03/'merge_active_batches.py',['--output',work/'merged_active'])
    run('evaluate_full360',p03/'evaluate_forward.py',['--freeze',freeze,'--active',active,'--initial-old-batch',raw/'formal_forward_v2_initial278','--initial-resumed-batch',raw/'formal_forward_v2_resume124','--active-old-batch',raw/'formal_forward_v2_active_new','--active-resumed-batch',raw/'formal_forward_v2_active_resume14','--initial-batch',raw/'formal_forward_v2_merged','--active-batch',raw/'formal_forward_v2_active_merged41','--initial-qc',r2/'05_independent_audit/formal_v2_initial_independent_qc.csv','--active-qc',r2/'05_independent_audit/formal_v2_active_independent_qc.csv','--initial-contacts',r2/'03_compute/method/formal_v2_initial_middle_contacts.csv','--active-contacts',r2/'03_compute/method/formal_v2_active_middle_contacts.csv','--output',work/'evaluation_replay'])
    for n in ('forward_arm_summary.csv','forward_attempt_level.csv','forward_charge_sensitivity.csv','forward_strata.csv'):
        if n=='forward_attempt_level.csv':scientific_table_same('R2_evaluation_'+n,ROOT/'data'/n,work/'evaluation_replay'/n,['independent_qc_evidence'])
        else:same('R2_evaluation_'+n,ROOT/'data'/n,work/'evaluation_replay'/n)
    run('small_example',ROOT/'run_example.py',['--output',work/'small_example']+([] if a.skip_small_training else ['--train-small']))
    plots=work/'plot_work';plots.mkdir(exist_ok=True);shutil.copytree(ROOT/'data',plots/'data',dirs_exist_ok=True);shutil.copytree(ROOT/'05_independent_audit',plots/'05_independent_audit',dirs_exist_ok=True);(plots/'figures').mkdir(exist_ok=True)
    shutil.copytree(ROOT/'plot_inputs',plots/'plot_inputs',dirs_exist_ok=True)
    shutil.copytree(ROOT/'plots',plots/'plots',dirs_exist_ok=True)
    run('plot_F8',plots/'plots/plot_forward_three_arms.py',extra={'PACKAGE_ROOT':str(plots)})
    run('plot_stages',plots/'plots/build_verified_stage_figures.py',extra={'PACKAGE_ROOT':str(plots)})
    run('plot_confirmation',plots/'plots/plot_from_verified_csv.py',['F4',plots/'data/confirmation_plot.csv',plots/'figures/confirmation_four_conditions'])
    chemical=work/'chemical_replay'
    if not chemical.exists():shutil.copytree(ROOT/'r3/library',chemical)
    run('chemical_domain9098',chemical/'scripts/prepare_frozen_domain.py',[freeze])
    run('chemical_alerts9098',chemical/'scripts/run_chemical_alerts.py')
    for n in ('chemical_alerts_9098.csv','alert_summary.csv','rule_frequencies.csv','figure_source.csv'):
        same('R3_chemical_'+n,ROOT/'r3/library/diagnostics'/n,chemical/'diagnostics'/n)
    anm=ROOT/'r3/anm';anm_config=work/'ANM_INPUTS.internal.json'
    anm_config.write_text(json.dumps({'H_AF_domain':inputs['human_aligned_pdb'],'H_AF_full_chainA':inputs['human_aligned_pdb'],'M_PUB_ALL_domain':inputs['mouse_aligned_pdb'],'M_PUB_ALL_full_chainA':inputs['mouse_aligned_pdb']},indent=2),encoding='utf-8')
    run('ANM_check',anm/'reproduce_anm.py',['--inputs',anm_config,'--out',work/'anm_replay','--check-only'])
    run('ANM_actual',anm/'reproduce_anm.py',['--inputs',anm_config,'--out',work/'anm_replay'])
    for n in ('all_eigenvalues.csv','low20_modes.csv','low20_relative_MSF.csv','boundary_sensitivity.csv','cutoff_sensitivity.csv','CA_node_index.csv','conformation_diagnostics.csv','conformation_displacements.csv'):
        same('R3_ANM_'+n,anm/n,work/'anm_replay/release'/n)
    run('ANM_figures',anm/'plot_anm_results.py',['--release',work/'anm_replay/release'])
    same('R3_ANM_mode_localization.csv',anm/'mode_localization.csv',work/'anm_replay/release/mode_localization.csv')
    run('structure_raw26',ROOT/'r3/structure/readback_qc.py',['--mouse',inputs['mouse_pdbqt'],'--human',inputs['human_pdbqt'],'--output',work/'structure_readback.csv'])
    s=json.loads((work/'structure_readback.json').read_text());assert s['jobs']==26 and s['protein_qc_passed']==25 and s['original_qc_agreement']==26
    refs=ROOT/'r3/reference_generation';refout=work/'reference_readback';refout.mkdir(exist_ok=True)
    for name,count in [('top10',9),('reference_panel',7),('denatonium',1)]:
        run('reference_'+name,refs/'audit_real_jobs_R3.py',[name,name+'_QC',refs/'raw_runs',Path(inputs['mouse_pdbqt']).parent],cwd=refs,extra={'TMC1_QC_OUTPUT_DIR':str(refout)})
        s=json.loads((refout/(name+'_QC_summary.json')).read_text());assert s['jobs_in_ledger']==s['independent_qc_pass']==count
    receipt={'status':'PASS','scope':'full existing replay: rankings, acquisition,360-attempt evaluation,small example,figures,9098 alerts,6 ANM networks/24 CA traversals,26 static pose readback,17 earlierR3 raw jobs (9 generated+7 paper references+1 denatonium)','actual_steps':len(receipts),'successful_unique_steps':len({r['id'] for r in receipts if r['exit_code']==0}),'failed_child_invocations_retained':sum(r['exit_code']!=0 for r in receipts),'full_ranking_run':not a.skip_full_ranking,'small_training_run':not a.skip_small_training,'original_scientific_components_changed':False,'new_docking_jobs':0,'new_MD_or_full_atom_dynamic_receptors':0,'byte_identity_checks':checks,'external_input_sha256':{k:sha(Path(v)) for k,v in inputs.items()},'logs_contain_local_paths':'independent work/internal only; publish sanitized receipts'}
    (work/'INTEGRATED_REPLAY_REPORT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('FULL EXISTING INTEGRATED REPLAY PASS',len(receipts),'steps',flush=True)
if __name__=='__main__':main()
