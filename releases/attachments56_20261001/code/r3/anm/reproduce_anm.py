"""Safe replay wrapper for a new output directory and external PDB inputs."""
import argparse,hashlib,json,os,subprocess,sys
from pathlib import Path
def main():
    root=Path(__file__).resolve().parent
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--check-only',action='store_true');a=p.parse_args()
    freeze=json.loads((root/'METHOD_FREEZE.json').read_text(encoding='utf-8'))
    protocol=root/'PROTOCOL_FROZEN.json';script=root/'run_anm.py';reference=root/'reference_code/ProDy_v2.6.1_anm.py'
    def sha(x):return hashlib.sha256(x.read_bytes()).hexdigest()
    for x,expected in [(protocol,freeze['protocol_sha256']),(script,freeze['run_script_sha256']),(reference,freeze['official_reference_sha256'])]:
        if sha(x)!=expected:raise ValueError('Frozen method changed: '+x.name)
    data=json.loads(a.inputs.read_text(encoding='utf-8'));data['prody_reference_source']=str(reference.resolve())
    for spec in freeze['inputs']:
        value=data.get(spec['id'],'')
        if not value or '<' in value:raise ValueError('Set external source PDB for '+spec['id'])
        x=Path(value).expanduser()
        if not x.is_absolute():x=a.inputs.resolve().parent/x
        if not x.is_file() or sha(x)!=spec['sha256']:raise ValueError('External source missing or wrong SHA: '+spec['id'])
        data[spec['id']]=str(x.resolve())
    out=a.out.resolve()
    if out==root or out.is_relative_to(root) or root.is_relative_to(out):raise ValueError('Replay output must be separate from the frozen package')
    if out.exists() and any(out.iterdir()):raise FileExistsError('Replay output must be new or empty')
    if a.check_only:print('PASS: method/source hashes checked; output directory independent; no calculation performed');return
    out.mkdir(parents=True,exist_ok=True);local=out/'EXTERNAL_INPUT_PATHS.internal.json';local.write_text(json.dumps(data,indent=2),encoding='utf-8')
    env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
    command=[sys.executable,str(script),'--protocol',str(protocol),'--inputs',str(local),'--out',str(out)]
    with (out/'ACTUAL_REPLAY_EXECUTION.log').open('w',encoding='utf-8') as h:r=subprocess.run(command,env=env,stdout=h,stderr=subprocess.STDOUT)
    print('Actual ANM replay exit code:',r.returncode);raise SystemExit(r.returncode)
if __name__=='__main__':main()
