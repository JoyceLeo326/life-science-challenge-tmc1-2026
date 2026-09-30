"""Rebuild frozen experimental worm receptor from official 7USX PDB."""
import argparse,subprocess,hashlib,json
from pathlib import Path
O=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--source-pdb',required=True,type=Path);ap.add_argument('--prepare-bin',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if sha(a.source_pdb)!='c0888d1eacef231d9ec4a53ba2aac582ab3cbc0176065b55009f7f7410a19caf':raise ValueError('Official7USX source hash mismatch')
a.output.mkdir(parents=True,exist_ok=True);src=a.output/'7USX_protein_native.pdb';src.write_text('\n'.join(l for l in a.source_pdb.read_text().splitlines() if l.startswith('ATOM  ') and l[76:78].strip() not in ['H','D'])+'\nEND\n',encoding='utf-8')
b=json.loads((O/'protocol/worm_X_box.json').read_text());prefix=a.output/'W_7USX_LOCAL_COMPLETE'
cmd=[str(a.prepare_bin),'--read_pdb',str(src),'-o',str(prefix),'-p','-j','--write_pdb',str(prefix)+'_H.pdb','--delete_bad_res_from_box_radius','8','--box_center',*[str(v) for v in b['center']],'--box_size',*[str(v) for v in b['size']]]
rc=subprocess.call(cmd)
if rc:raise SystemExit(rc)
actual=sha(prefix.with_suffix('.pdbqt'));expected='907ef675c3f92bbb3efbec7234cbbbf8a3fb1bd412f4191ce574e238ae9b6651';print(json.dumps({'actual_sha256':actual,'expected_sha256':expected,'match':actual==expected}));raise SystemExit(0 if actual==expected else 2)
