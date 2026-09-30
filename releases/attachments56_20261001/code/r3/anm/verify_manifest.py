"""Verify the explicit anonymous NMA publication file manifest."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent
j=json.loads((root/'MANIFEST.json').read_text(encoding='utf-8'))
for row in j['files']:
    p=root/row['path']
    if not p.is_file() or p.stat().st_size!=row['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Manifest mismatch: '+row['path'])
print('MANIFEST PASS:',len(j['files']),'/',len(j['files']))
