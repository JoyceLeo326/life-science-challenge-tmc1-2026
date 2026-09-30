"""Read frozen pose files only; no docking. Supply external mouse/human receptors."""
import argparse,csv,json,hashlib,re,math
from pathlib import Path
import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree
ap=argparse.ArgumentParser();ap.add_argument('--mouse',required=True,type=Path);ap.add_argument('--human',required=True,type=Path);ap.add_argument('--worm',type=Path);ap.add_argument('--output',type=Path);a=ap.parse_args()
O=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(p):return list(csv.DictReader(p.open(encoding='utf-8-sig')))
def atoms(p):
 out=[]
 for l in p.read_text().splitlines():
  if l.startswith('ENDMDL'):break
  if l.startswith(('ATOM  ','HETATM')) and l.split()[-1] not in ['H','HD','HS'] and not re.fullmatch(r'G\d+',l.split()[-1]):out.append((l[21],int(l[22:26]),np.array([float(l[t:t+8]) for t in [30,38,46]])))
 return out
def stable(m,ns):
 m=Chem.RemoveHs(m)
 for n in ns:m.GetAtomWithIdx(n).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
 Chem.AssignStereochemistry(m,cleanIt=True,force=True);return Chem.MolToSmiles(m,isomericSmiles=True)
receptors={'M_PUB_ALL':a.mouse,'H_AF':a.human,'W_7USX_LOCAL_COMPLETE':a.worm or O/'protocol/W_7USX_LOCAL_COMPLETE.pdbqt'}
sites={'M_PUB_ALL':{408,411,412,443,444,447,448,528,531,532,536,579,580,582,601},'H_AF':{414,417,418,449,450,453,454,531,534,535,539,582,583,585,604},'W_7USX_LOCAL_COMPLETE':{400,403,404,435,436,439,440,683,686,687,691,734,735,737,756}}
cache={};out=[]
for j in rows(O/'all26_jobs_anonymous.csv'):
 bad=[];rec=receptors[j['receptor_id']]
 if sha(rec)!=j['receptor_sha256']:bad.append('receptor_hash')
 for field,hsh in [('ligand_pdbqt','ligand_sha256'),('sdf_file','sdf_sha256'),('pose_file','pose_sha256'),('exported_sdf_file','exported_sdf_sha256'),('log_file','log_sha256'),('box_file','box_sha256')]:
  if sha(O/j[field])!=j[hsh]:bad.append(field+'_hash')
 if str(rec) not in cache:
  ats=atoms(rec);cache[str(rec)]=(cKDTree(np.array([z[2] for z in ats])),np.array([z[2] for z in ats if z[0]=='A' and z[1] in sites[j['receptor_id']]]))
 tr,mid=cache[str(rec)];xyz=np.array([z[2] for z in atoms(O/j['pose_file'])]);ns=json.loads(j.get('exchangeable_protonated_N_indices') or '[]');rule=Chem.MolFromSmiles(j['rule_smiles'])
 for field in ['sdf_file','exported_sdf_file']:
  mol=next(Chem.SDMolSupplier(str(O/j[field]),removeHs=False))
  if stable(mol,ns)!=stable(rule,ns):bad.append(field+'_stable_graph')
 if len(xyz)!=rule.GetNumHeavyAtoms():bad.append('heavy_count')
 text=(O/j['pose_file']).read_text();score=float(re.search(r'REMARK VINA RESULT:\s*([-\d.]+)',text).group(1));log=(O/j['log_file']).read_text(encoding='latin-1')
 tab=float(next(l.split()[1] for l in log.splitlines() if re.match(r'^\s*1\s+-?\d',l)))
 if abs(score-float(j['score_kcal_mol']))>1e-6 or abs(score-tab)>.0051:bad.append('score_agreement')
 mind=float(tr.query(xyz)[0].min());margin=float(np.min(np.array(json.loads(j['box_size_A']))/2-np.abs(xyz-np.array(json.loads(j['box_center_A'])))))
 if mind<1.5:bad.append('protein_clash')
 if margin<-.1:bad.append('box_margin')
 md=float(np.linalg.norm(xyz[:,None,:]-mid[None,:,:],axis=2).min())
 expected=str(j['independent_qc_pass']).lower()=='true';actual=not bad
 out.append({'run_id':j['run_id'],'receptor_id':j['receptor_id'],'readback_protein_qc_pass':actual,'original_independent_protein_qc_pass':expected,'readback_agrees_with_original_qc':actual==expected,'failures':';'.join(bad),'score_kcal_mol':score,'min_protein_A':mind,'min_box_margin_A':margin,'closest_middle_A':md,'middle_contact4A':md<4})
dest=a.output or O/'public_readback_qc.csv'
with dest.open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,list(out[0]));w.writeheader();w.writerows(out)
summary={'jobs':len(out),'protein_qc_passed':sum(r['readback_protein_qc_pass'] for r in out),'protein_qc_failed':sum(not r['readback_protein_qc_pass'] for r in out),'original_qc_agreement':sum(r['readback_agrees_with_original_qc'] for r in out),'unexpected_disagreements':sum(not r['readback_agrees_with_original_qc'] for r in out),'scope':'frozen protein-only poses; original worm HET clashes separately preserved; known worm indinavir box failure retained; no score threshold'}
dest.with_suffix('.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary));raise SystemExit(0 if summary['unexpected_disagreements']==0 else 2)
