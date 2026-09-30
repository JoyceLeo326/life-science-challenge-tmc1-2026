"""Actual new docking replay. Do not overwrite frozen evidence."""
import argparse,os,subprocess,sys,hashlib
from pathlib import Path
O=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--species',choices=['worm','human'],required=True);ap.add_argument('--vina-bin',required=True,type=Path);ap.add_argument('--export-bin',required=True,type=Path);ap.add_argument('--receptor-dir',type=Path);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--limit',type=int,default=0);a=ap.parse_args()
if a.output.resolve()==O or O in a.output.resolve().parents:raise ValueError('Replay output must be outside frozen release folder')
worm=a.species=='worm';rec='W_7USX_LOCAL_COMPLETE' if worm else 'H_AF';rd=a.receptor_dir or O/'protocol'
if not worm and a.receptor_dir is None:raise ValueError('Human H_AF author model is not redistributed; supply exact authorized external receptor directory')
os.environ['TMC1_R1_ROOT']=str(O);os.environ['TMC1_WORK_ROOT']=str(a.output.resolve());os.environ['PYTHONDONTWRITEBYTECODE']='1'
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
cmd=[sys.executable,'-B',str(O/'protocol'/('dock_worm_batch.py' if worm else 'dock_human_reference2.py')),'--input',str(O/('frozen_worm_input11.csv' if worm else 'frozen_human_reference2.csv')),'--batch','replay_'+a.species,'--scratch',str(a.output),'--receptor',rec,'--receptor-dir',str(rd),'--vina-bin',str(a.vina_bin),'--export-bin',str(a.export_bin),'--box',str(O/'protocol'/('worm_X_vina_box.txt' if worm else 'human_vina_box.txt')),'--box-meta',str(O/'protocol'/('worm_X_box.json' if worm else 'human_common_box.json')),'--workers','2','--seed','20260927','--exhaustiveness','8','--timeout','1800']
if a.limit:cmd+=['--limit',str(a.limit)]
raise SystemExit(subprocess.call(cmd))
