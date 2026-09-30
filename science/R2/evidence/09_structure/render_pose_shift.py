"""Render independently QC-passed seed29 parent and design SDF poses in one frame."""
from pathlib import Path
import os
import hashlib
import json
import time
import pymol2

HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = Path(os.environ.get('PACKAGE_ROOT', Path(__file__).resolve().parents[1]))
ROOT = PACKAGE_ROOT
REC = Path(r'<PACKAGE_ROOT>')
PARENT = ROOT / r'03_算法与计算\raw_runs\design_parent_pair\poses\LIB_CHEMBL3894860__M_PUB_ALL__e8__s20260929.sdf'
DESIGN = ROOT / r'03_算法与计算\raw_runs\design_fixed_remaining5\poses\DES_TWTKMJWAXOPXSX__M_PUB_ALL__e8__s20260929.sdf'
OUT = Path(r'<PACKAGE_ROOT>\Documents\Codex\TMC1_R2_Scratch_20260929\09_visuals\panel_D_pose_shift.png')
MIDDLE = [408,411,412,443,444,447,448,528,531,532,536,579,580,582,601]

def main():
    t=time.time()
    # Vina SDF contains nine ranked poses; restrict the display to rank 1.
    def first_pose(path, label):
        block=path.read_text(encoding='utf-8',errors='replace').split('$$$$',1)[0]+'$$$$\n'
        target=OUT.parent/(label+'_seed29_first_pose.sdf')
        target.write_text(block,encoding='utf-8')
        return target
    parent_first=first_pose(PARENT,'parent')
    design_first=first_pose(DESIGN,'design')
    with pymol2.PyMOL() as pm:
        c=pm.cmd
        for k,v in {'max_threads':1,'orthoscopic':1,'antialias':1,'ray_shadows':0,
                    'depth_cue':0,'specular':0.16,'ambient':0.55,'direct':0.70,
                    'stick_radius':0.29,'sphere_scale':0.37}.items(): c.set(k,v)
        c.bg_color('white')
        c.load(str(REC),'rec')
        c.load(str(parent_first),'parent')
        c.load(str(design_first),'design')
        atom_counts={'parent':c.count_atoms('parent and not hydro'), 'design':c.count_atoms('design and not hydro')}
        c.dss('rec')
        c.hide('everything','all')
        c.select('focus','byres (rec and chain A within 12 of (parent or design))')
        c.show('cartoon','focus')
        c.color('gray80','focus')
        c.set('cartoon_transparency',0.28,'focus')
        c.select('middle','rec and chain A and resi '+'+'.join(map(str,MIDDLE)))
        c.show('spheres','middle and name CA')
        c.color('forest','middle and name CA')
        c.set('sphere_scale',0.49,'middle and name CA')
        for obj,rgb in [('parent',(0.02,0.47,0.65)),('design',(0.87,0.32,0.20))]:
            c.set_color(obj+'_carbon',rgb)
            c.show('sticks',obj)
            c.color(obj+'_carbon',obj+' and elem C')
            c.color('atomic',obj+' and not elem C')
            c.set('stick_radius',0.34,obj)
        c.orient('focus or parent or design')
        c.turn('y',15)
        c.turn('x',10)
        c.zoom('focus or parent or design',buffer=12)
        view=list(c.get_view())
        c.png(str(OUT),width=1500,height=1100,dpi=300,ray=1,quiet=1)
        version=c.get_version()[0]
    meta={'source_coordinate_file':str(REC),'source_sha256':hashlib.sha256(REC.read_bytes()).hexdigest(),
          'parent_pose':str(PARENT),'parent_pose_sha256':hashlib.sha256(PARENT.read_bytes()).hexdigest(),
          'parent_rank1_sha256':hashlib.sha256(parent_first.read_bytes()).hexdigest(),
          'design_pose':str(DESIGN),'design_pose_sha256':hashlib.sha256(DESIGN.read_bytes()).hexdigest(),
          'design_rank1_sha256':hashlib.sha256(design_first.read_bytes()).hexdigest(),
          'seed':20260929,'receptor':'M_PUB_ALL','box':'table6_middle_15_common_union',
          'middle_residues_mouse_chain_A':MIDDLE,'heavy_atom_counts_loaded':atom_counts,
          'camera_view':view,'renderer':'PyMOL '+str(version),'elapsed_seconds':round(time.time()-t,2)}
    (HERE/'pose_shift_render_manifest.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'output':str(OUT),'atoms':atom_counts,'seconds':meta['elapsed_seconds']},ensure_ascii=False))

if __name__=='__main__': main()
