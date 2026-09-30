"""Real small design/optimization/screening example using package inputs.

This separate tutorial does not alter the frozen six-candidate result list.
"""
import argparse,csv,json,os,subprocess,sys,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
    p=argparse.ArgumentParser();p.add_argument('--bundle',type=Path,default=ROOT/'reproduction_bundle');p.add_argument('--output',type=Path,default=ROOT/'example_results');p.add_argument('--train-small',action='store_true');a=p.parse_args()
    bundle=a.bundle.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    evidence=bundle/'science/R2/evidence';method=out/'_method'
    shutil.copytree(evidence/'04_design/method_results',method,dirs_exist_ok=True)
    env=os.environ.copy();env.update(TMC1_R1_ROOT=str(bundle/'science/R1'),TMC1_R2_ROOT=str(evidence),TMC1_WORK_ROOT=str(out),TMC1_RULE_OUTPUT=str(out/'rule_design'))
    (out/'rule_design').mkdir(exist_ok=True)
    commands=[[sys.executable,str(method/'run_design.py')]]
    if a.train_small:
        for script,extra,runid in [('train_smiles_gru.py',['--epochs','1','--samples','32','--batch-size','16'],'example_smiles'),('train_selfies_gru.py',['--base-epochs','1','--finetune-epochs','0','--samples','32'],'example_selfies')]:commands.append([sys.executable,str(method/script),'--library',str(ROOT/'data/tutorial_library2048.csv'),'--smiles-column','standardized_isomeric_smiles','--run-id',runid]+extra)
    receipts=[]
    for i,cmd in enumerate(commands):
        log=out/('step'+str(i)+'.log')
        with log.open('wb') as f:r=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
        receipts.append({'entry':Path(cmd[1]).name,'exit_code':r.returncode,'log':log.name})
        if r.returncode:raise RuntimeError('Failed '+Path(cmd[1]).name+'; inspect '+log.name)
    # Optimize via the original frozen property/alert constraints, then apply a
    # transparent illustrative neighborhood prior; no activity model claim.
    table=list(csv.DictReader((out/'rule_design/parent_product_comparison_data.csv').open(encoding='utf-8-sig')))
    fixed={r['compound_id'] for r in csv.DictReader((out/'rule_design/fixed_docking_list.csv').open(encoding='utf-8-sig'))}
    selected=[r for r in table if r['role']=='enumerated_analogue' and r['compound_id'] in fixed]
    assert len(selected)==7
    selected.sort(key=lambda r:(float(r['docking_neighborhood_prior_kcal_mol']),r['compound_id']))
    for i,r in enumerate(selected,1):r.update(example_rank=i,score_type='uncalibrated parent-neighborhood prior; not affinity',mode='separate reproducible rule-enumeration example')
    with (out/'example_screening.csv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,selected[0].keys());w.writeheader();w.writerows(selected)
    (out/'EXAMPLE_RECEIPT.json').write_text(json.dumps({'actual_steps':receipts,'designs':len(selected),'formal_results_modified':False},indent=2),encoding='utf-8')
    print('Example completed: actual rule design -> property/alert selection -> illustrative screening',len(selected),'designs')
if __name__=='__main__':main()
