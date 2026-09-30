"""Register actual plotted files and sources without inventing visual approval.

Adapted for the final relative directory layout. A human/visual approval is a
separate state; this default command verifies file integrity and records it.
"""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    records=[]
    for png in sorted((ROOT/'figures').glob('*.png')):
        pdf=png.with_suffix('.pdf');svg=png.with_suffix('.svg')
        if not pdf.is_file() or not svg.is_file():raise FileNotFoundError('Incomplete figure outputs: '+png.stem)
        records.append({'figure':png.stem,'files':[{'path':f.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in [png,pdf,svg]],'visual_checked':False,'status':'rendered and file-integrity verified; visual approval remains separate'})
    if not records:raise ValueError('No rendered figures')
    (ROOT/'figures/FIGURE_INDEX.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Registered',len(records),'real rendered figures')
if __name__=='__main__':main()
