"""Run a frozen R3 batch with user-provided Vina1.2.7 and Meeko export executable."""
from pathlib import Path
import argparse,sys,os,subprocess,shutil,json,hashlib
p=argparse.ArgumentParser()
p.add_argument('--batch',choices=['top10','reference_panel','denatonium'],default='top10')
p.add_argument('--vina',required=True,type=Path)
p.add_argument('--export',required=True,type=Path,help='Meeko mk_export executable, e.g. venv/Scripts/mk_export.exe')
p.add_argument('--receptor-dir',required=True,type=Path,help='Externally supplied M_PUB_ALL.pdbqt; third-party receptor is not redistributed')
p.add_argument('--limit',type=int,default=0,help='0 full frozen batch; 1 prepares/docks only the first selected row for a light entry check')
p.add_argument('--out',type=Path,default=Path('R3_replay_new_runs'))
a=p.parse_args();root=Path(__file__).resolve().parent
receptor=a.receptor_dir.resolve()/'M_PUB_ALL.pdbqt'
expected='c024cc9efaf910ac11642cea4a4c701a0a88a335d7deb0d040320516d5a2eda1'
if not receptor.is_file():p.error('External M_PUB_ALL.pdbqt is required; it is not included in this release.')
if hashlib.sha256(receptor.read_bytes()).hexdigest()!=expected:p.error('External receptor SHA256 differs from the frozen protocol.')
if a.limit<0:p.error('--limit must be nonnegative')
out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
name={'top10':'frozen_top10.csv','reference_panel':'reference_identity_inputs.csv','denatonium':'denatonium_input.csv'}[a.batch]
env=os.environ.copy();env.update(TMC1_R1_ROOT=str(root),TMC1_WORK_ROOT=str(out),OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',TMC1_QC_OUTPUT_DIR=str(out/'qc'))
cmd=[sys.executable,str(root/'protocol_inputs/dock_batch.py'),'--input',str(root/name),'--batch',a.batch,'--scratch',str(out/'raw_runs'),'--receptor','M_PUB_ALL','--receptor-dir',str(a.receptor_dir.resolve()),'--vina-bin',str(a.vina.resolve()),'--export-bin',str(a.export.resolve()),'--box',str(root/'protocol_inputs/vina_box.txt'),'--box-meta',str(root/'protocol_inputs/common_box.json'),'--seed','20260927','--exhaustiveness','8','--workers','2','--timeout','600','--limit',str(a.limit)]
(out/'replay_command.json').write_text(json.dumps(cmd,indent=2),encoding='utf-8')
with (out/'replay_execution.log').open('w',encoding='utf-8') as log:subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
(out/'qc').mkdir(exist_ok=True)
subprocess.run([sys.executable,str(root/'audit_real_jobs_R3.py'),a.batch,'REPLAY_'+a.batch+'_QC',str(out/'raw_runs')],env=env,check=True)
