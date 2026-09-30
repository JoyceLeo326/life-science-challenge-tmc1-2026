"""Freeze aligned receptors and one common literature-defined search volume.

Rigid coordinate transforms only. No model completion, energy minimization or
humanization is performed. The author's complex is a public predicted/MD model.
"""
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.SeqUtils import seq1
from analyze_structure import FULL_SITE_TEXT
from common import ROOT, write_json, write_csv, sha256, run_logged

OUT=ROOT/'screen_v2/inputs/receptors'

def atoms(path):
    result=[]
    for line in path.read_text().splitlines():
        if not line.startswith(('ATOM  ','HETATM')): continue
        elem=line[76:78].strip() or line[12:16].strip()[0]
        if elem in ('H','D'):continue
        result.append(dict(line=line,name=line[12:16].strip(),chain=line[21],n=int(line[22:26]),
                           res=line[17:20].strip(),element=elem,xyz=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])],float)))
    return result

def ca_map(a,chain='A'):
    return {x['n']:x for x in a if x['name']=='CA' and x['chain']==chain}

def align(source,target,tm):
    sc=ca_map(source);tc=ca_map(target)
    nums=list(sc);hnums=list(tc)
    ss=''.join(seq1(x['res']) for x in sc.values());hs=''.join(seq1(x['res']) for x in tc.values())
    aln=Align.PairwiseAligner(mode='global',substitution_matrix=substitution_matrices.load('BLOSUM62'),open_gap_score=-10,extend_gap_score=-0.5).align(ss,hs)[0]
    pairs=[]
    for (a,b),(c,d) in zip(*aln.aligned):
        pairs.extend((nums[i],hnums[j]) for i,j in zip(range(a,b),range(c,d)))
    fit=[(s,h) for s,h in pairs if h in tm and sc[s]['res']==tc[h]['res']]
    x=np.array([sc[s]['xyz'] for s,h in fit]);y=np.array([tc[h]['xyz'] for s,h in fit]);xc=x.mean(0);yc=y.mean(0)
    u,_,vt=np.linalg.svd((x-xc).T@(y-yc));rot=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
    rmsd=float(np.sqrt(np.mean(np.sum(((x-xc)@rot+yc-y)**2,axis=1))))
    return pairs,rot,xc,yc,rmsd,fit

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    uniprot=json.loads((ROOT/'data/raw/uniprot_Q8TDI8.json').read_text())
    tm={n for f in uniprot['features'] if f['type']=='Transmembrane' for n in range(f['location']['start']['value'],f['location']['end']['value']+1)}
    sources=[('H_AF','structure_bundle/raw/AFDB/Q8TDI8/AF-Q8TDI8-F1.pdb','human','AFDB v6 canonical Q8TDI8 monomer'),
             ('M_AF','structure_bundle/raw/AFDB/Q8R4P5/AF-Q8R4P5-F1.pdb','mouse','AFDB v6 canonical Q8R4P5 monomer'),
             ('M_PUB_ALL','data/raw/public_MET_complex.pdb','mouse','Published authors AF2/MD six-protein coordinate model; no lipids/ions/water in supplied file')]
    human=atoms(ROOT/sources[0][1]);records=[];mapping_rows=[];all_site_xyz=[];preflight=[]
    mouse_to_human={int(r['mouse_resnum']):int(r['human_resnum']) for r in __import__('csv').DictReader((ROOT/'results/mouse_human_residue_map.csv').open(encoding='utf-8-sig'))}
    middle_human={mouse_to_human[int(t[1:])] for t in FULL_SITE_TEXT['middle'].split()}
    for rid,rel,species,description in sources:
        raw=atoms(ROOT/rel)
        pairs,rot,xc,yc,rmsd,fit=align(raw,human,tm)
        sn_to_hn=dict(pairs)
        lines=[];near_site=[];ca=ca_map(raw)
        for a in raw:
            xyz=(a['xyz']-xc)@rot+yc
            lines.append(a['line'][:30]+''.join(f'{v:8.3f}' for v in xyz)+a['line'][54:])
            if a['chain']=='A' and sn_to_hn.get(a['n']) in middle_human:
                near_site.append(xyz)
                all_site_xyz.append(xyz)
        mapped_middle={sn_to_hn[n] for n in ca if sn_to_hn.get(n) in middle_human}
        if mapped_middle!=middle_human:raise ValueError(f'{rid} missing middle-site residues: {middle_human-mapped_middle}')
        aligned=OUT/f'{rid}_aligned_heavy.pdb';aligned.write_text('\n'.join(lines)+'\nEND\n',encoding='utf-8')
        for s,h in pairs:
            mapping_rows.append({'receptor_id':rid,'source_chain':'A','source_resnum':s,'source_aa':seq1(ca[s]['res']),
                'human_resnum':h,'human_aa':seq1(ca_map(human)[h]['res']),'conserved':ca[s]['res']==ca_map(human)[h]['res'],
                'in_tm_fit':(s,h) in fit,'in_table6_middle':h in middle_human})
        sg=[a for a in raw if a['name']=='SG'];disulfides=[]
        for i,a in enumerate(sg):
            for b in sg[i+1:]:
                distance=float(np.linalg.norm(a['xyz']-b['xyz']))
                if distance<2.5:disulfides.append({'residue1':f"{a['chain']}:{a['n']}",'residue2':f"{b['chain']}:{b['n']}",'sg_distance_A':distance,'evidence':'model geometry only, not experimental disulfide validation'})
        preflight.append({'receptor_id':rid,'source':rel,'source_sha256':sha256(ROOT/rel),'source_heavy_atoms':len(raw),'source_chains':dict(__import__('collections').Counter(a['chain'] for a in raw)),
            'chainA_ca_count':len(ca),'aligned_pairs':len(pairs),'conserved_tm_fit_pairs':len(fit),'tm_ca_rmsd_A':rmsd,
            'rotation':rot.tolist(),'source_centroid':xc.tolist(),'target_centroid':yc.tolist(),'geometric_disulfide_pairs':disulfides,
            'middle_residues_present':len(mapped_middle),'alignment_effect':'Rigid transform; source internal geometry unchanged apart from PDB 0.001 A rounding'})
        prefix=f'screen_v2/inputs/receptors/{rid}'
        run_logged([Path(sys.executable).parent/'mk_prepare_receptor.exe','--read_pdb',str(aligned.relative_to(ROOT)),'-o',prefix,'-p','-j','--write_pdb',prefix+'_H.pdb'],f'v2_prepare_{rid}')
        records.append({'receptor_id':rid,'species':species,'description':description,'source':rel,'source_sha256':sha256(ROOT/rel),
                        'aligned_pdb':str(aligned.relative_to(ROOT)),'pdbqt':prefix+'.pdbqt','pdbqt_sha256':sha256(ROOT/(prefix+'.pdbqt')),
                        'interpretation':'model-condition sensitivity; no experimentally assigned open/closed state'})
    # Controlled context ablation is extracted from the same already-parameterized complex.
    # This prevents a second preparation from changing atom types or coordinate order in TMC1.
    full=(OUT/'M_PUB_ALL.pdbqt').read_text().splitlines()
    core=[l for l in full if l.startswith(('ATOM  ','HETATM')) and l[21]=='A']
    (OUT/'M_PUB_CORE.pdbqt').write_text('\n'.join(core)+'\n',encoding='utf-8')
    if not core:raise ValueError('No chain A after receptor preparation')
    core_heavy=[l for l in (OUT/'M_PUB_ALL_aligned_heavy.pdb').read_text().splitlines() if l.startswith('ATOM') and l[21]=='A']
    (OUT/'M_PUB_CORE_aligned_heavy.pdb').write_text('\n'.join(core_heavy)+'\nEND\n',encoding='utf-8')
    records.append({'receptor_id':'M_PUB_CORE','species':'mouse','description':'Chain A subset of identical parameterized M_PUB_ALL; controlled protein-context ablation',
        'source':'screen_v2/inputs/receptors/M_PUB_ALL.pdbqt','source_sha256':sha256(OUT/'M_PUB_ALL.pdbqt'),
        'aligned_pdb':'screen_v2/inputs/receptors/M_PUB_CORE_aligned_heavy.pdb','pdbqt':'screen_v2/inputs/receptors/M_PUB_CORE.pdbqt',
        'pdbqt_sha256':sha256(OUT/'M_PUB_CORE.pdbqt'),'interpretation':'Removing other protein chains only; not isolating lipids, membrane or gating'})
    xyz=np.array(all_site_xyz);lo=xyz.min(0)-5;hi=xyz.max(0)+5
    box={'site_id':'table6_middle_15_common_union','center':((lo+hi)/2).round(4).tolist(),'size':(hi-lo).round(4).tolist(),
         'volume_A3':float(np.prod(hi-lo)),'human_site_residues':sorted(middle_human),'mouse_site_residues':FULL_SITE_TEXT['middle'].split(),
         'rule':'Union of all heavy atoms of the 15 Supplementary Table 6 middle-site mapped residues in H_AF, M_AF, M_PUB_ALL chain A, after conserved-TM rigid alignment; 5 A padding per face; SAME center/size for every receptor and ligand',
         'source_conflict':'Paper body calls mouse M407 middle, Supplementary Table 6 calls it Top. Project pre-specifies Table 6 region labels; M407 is not in this 15-residue middle set.',
         'interpretation':'Literature-informed regional hypothesis, not a validated human binding pocket or gate state'}
    write_json(ROOT/'screen_v2/config/common_box.json',box)
    (ROOT/'screen_v2/config/vina_box.txt').write_text('\n'.join(f'{p}_{ax} = {v:.4f}' for p,vals in [('center',box['center']),('size',box['size'])] for ax,v in zip('xyz',vals))+'\n')
    for record in records:
        a=atoms(ROOT/record['aligned_pdb']);inside=[x for x in a if np.all(np.abs(x['xyz']-box['center'])<=np.array(box['size'])/2+8)]
        record['heavy_atoms_within_box_plus8A_by_chain']=dict(__import__('collections').Counter(x['chain'] for x in inside))
    write_json(OUT/'receptors.json',records);write_csv(OUT/'residue_mapping.csv',mapping_rows);write_json(OUT/'structural_preflight.json',preflight)
    print(json.dumps({'receptor_count':len(records),'box':box,'context_counts':[(r['receptor_id'],r['heavy_atoms_within_box_plus8A_by_chain']) for r in records]},indent=2),flush=True)

if __name__=='__main__':main()
