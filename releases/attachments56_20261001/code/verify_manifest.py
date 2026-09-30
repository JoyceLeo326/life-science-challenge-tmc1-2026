"""Verify every file listed in the final publication manifest."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
    j=json.loads((ROOT/'MANIFEST.json').read_text(encoding='utf-8'));errors=[]
    for r in j['files']:
        p=ROOT/r['path']
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:errors.append(r['path'])
    if errors:raise ValueError('Manifest failures: '+repr(errors))
    print('MANIFEST PASS:',len(j['files']),'/',len(j['files']))
if __name__=='__main__':main()
