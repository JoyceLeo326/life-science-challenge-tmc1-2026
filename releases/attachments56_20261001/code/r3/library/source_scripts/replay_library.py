"""Offline reconstruction into a new independent directory, never the public input tree."""
import argparse,os,shutil,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--python',default=sys.executable);p.add_argument('--check-only',action='store_true');a=p.parse_args()
base=Path(__file__).resolve().parent.parent;out=a.out.resolve()
if out==base or base in out.parents:raise SystemExit('Output must be outside this release')
if out.exists():raise SystemExit('Output must be a new directory')
raw=base/'provenance/raw/expanded';legacy=base/'provenance/raw/legacy_approved';history=base/'provenance/historical'
if len(list(raw.glob('offset_*.json')))!=20 or len(list(legacy.glob('molecules_offset_*.json')))!=4:raise SystemExit('Expected 20+4 raw pages')
if a.check_only:print('24 raw page files present; no scientific processing performed');raise SystemExit(0)
out.mkdir();work=out/'source_scripts';work.mkdir()
for name in ['build_library.py','make_overlap.py','verify_and_profile.py']:shutil.copyfile(base/'source_scripts'/name,work/name)
env=dict(os.environ,TMC1_LIBRARY_WORK=str(work),TMC1_R2_RAW_CACHE=str(raw),TMC1_R1_APPROVED_CACHE=str(legacy),TMC1_R1_ROOT=str(history),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
for cmd in [['build_library.py','prepare'],['make_overlap.py'],['verify_and_profile.py']]:subprocess.run([a.python,str(work/cmd[0]),*cmd[1:]],check=True,env=env,cwd=work)
print('New independent outputs:',work/'processed')
