from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFEST.json').read_text(encoding='utf-8'))
for row in manifest['files']:
 p=(root/row['path']).resolve()
 if not p.is_relative_to(root):raise ValueError('Manifest path escaped release')
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 if p.stat().st_size!=row['bytes'] or h.hexdigest()!=row['sha256']:raise ValueError('Mismatch: '+row['path'])
print('MANIFEST PASS:',len(manifest['files']),'/',len(manifest['files']))
