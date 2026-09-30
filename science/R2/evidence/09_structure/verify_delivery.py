"""Verify the deliverable figures and record hashes/dimensions for integration."""
from pathlib import Path
from PIL import Image
import csv
import hashlib
import json
import xml.etree.ElementTree as ET

here=Path(__file__).resolve().parent
stems=['图09_六链结构与双区域','图09_亲本设计物区域迁移']
files={}
for stem, expected in zip(stems,[(4800,2520),(4500,2310)]):
    png=here/(stem+'.png'); svg=here/(stem+'.svg')
    with Image.open(png) as im:
        assert im.size==expected and im.mode in ('RGB','RGBA')
    root=ET.parse(svg).getroot()
    ns={'s':'http://www.w3.org/2000/svg'}
    text_nodes=root.findall('.//s:text',ns)
    image_nodes=root.findall('.//s:image',ns)
    assert len(text_nodes)>=10 and len(image_nodes)>=1
    files[stem]={'png':{'path':str(png),'sha256':hashlib.sha256(png.read_bytes()).hexdigest(),
                        'dimensions_px':expected},
                 'svg':{'path':str(svg),'sha256':hashlib.sha256(svg.read_bytes()).hexdigest(),
                        'editable_text_nodes':len(text_nodes),'embedded_structure_images':len(image_nodes)}}
with (here/'图09_区域迁移数据.csv').open(encoding='utf-8-sig',newline='') as f:
    rows=list(csv.DictReader(f))
assert len(rows)==3 and all(r['design_independent_qc_pass']=='True' and r['parent_independent_qc_pass']=='True' for r in rows)
files['source_table']={'rows':len(rows),'sha256':hashlib.sha256((here/'图09_区域迁移数据.csv').read_bytes()).hexdigest()}
(here/'delivery_manifest.json').write_text(json.dumps(files,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'figures':len(stems),'three_seed_rows':len(rows),'checks':'pass'},ensure_ascii=False))
