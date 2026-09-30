"""Assemble existing frozen ligand poses with an authorized external receptor.

No docking, minimization, prediction, charge reassignment or coordinate transform.
Output PDB preserves receptor heavy coordinates and rounds SDF coordinates to
the PDB 0.001 Angstrom field precision. Original SDF retains exact bond orders.
"""
import argparse,csv,json,hashlib,shutil
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
from rdkit import Chem,rdBase

RECEPTOR_PDB_SHA='a0c2c7b07eb3a4939397ffb104d0d8b8bdc5ef4ab9b7af732deef03944385baa'
RECEPTOR_PDBQT_SHA='c024cc9efaf910ac11642cea4a4c701a0a88a335d7deb0d040320516d5a2eda1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def csvrows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def writecsv(p,rows):
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
def xyz(line):return tuple(float(line[i:i+8]) for i in (30,38,46))
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--receptor-pdb',required=True,type=Path);ap.add_argument('--receptor-pdbqt',required=True,type=Path);ap.add_argument('--candidates',type=Path,default=Path(__file__).resolve().parent/'results.csv');ap.add_argument('--ligands',type=Path,default=Path(__file__).resolve().parent/'structures');ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
    assert sha(a.receptor_pdb)==RECEPTOR_PDB_SHA,'External aligned heavy receptor SHA mismatch'
    assert sha(a.receptor_pdbqt)==RECEPTOR_PDBQT_SHA,'External frozen M_PUB_ALL receptor SHA mismatch'
    if (a.output/'ASSEMBLY_PROVENANCE.json').exists():raise FileExistsError('Completed assembly cannot be overwritten')
    receptor=[l.rstrip() for l in a.receptor_pdb.read_text().splitlines() if l.startswith(('ATOM  ','HETATM'))]
    assert len(receptor)==15070 and {l[21] for l in receptor}==set('ABCDEF')
    assert not any(l[76:78].strip()=='H' for l in receptor)
    qheavy=[l for l in a.receptor_pdbqt.read_text().splitlines() if l.startswith(('ATOM  ','HETATM')) and l.split()[-1] not in ('H','HD','HS')]
    assert len(qheavy)==len(receptor) and Counter(xyz(l) for l in receptor)==Counter(xyz(l) for l in qheavy),'Receptor PDB/PDBQT coordinate frame must match exactly'
    candidates=csvrows(a.candidates);assert len(candidates)==len({r['candidate_id'] for r in candidates})==6
    a.output.mkdir(parents=True,exist_ok=True);(a.output/'ligands').mkdir(exist_ok=True);(a.output/'complexes').mkdir(exist_ok=True)
    shutil.copyfile(a.candidates,a.output/'frozen_candidates.csv')
    index=[];atommap=[];bonds=[]
    for candidate in candidates:
        cid=candidate['candidate_id'];source=a.ligands/(cid+'.sdf')
        supplier=Chem.SDMolSupplier(str(source),removeHs=False);mol=supplier[0];assert mol is not None and len(supplier)>=1
        conf=mol.GetConformer();assert conf.Is3D()
        title=mol.GetProp('_Name');assert 'M_PUB_ALL' in title and 'e8' in title and 's20260927' in title
        sourcecopy=a.output/'ligands'/source.name;shutil.copyfile(source,sourcecopy)
        lines=['HEADER    COMPUTATIONAL POSE ASSEMBLY',
               'TITLE     '+cid+' WITH FROZEN M_PUB_ALL RECEPTOR',
               'REMARK 900 EXISTING FIRST VINA POSE; NO NEW DOCKING OR MINIMIZATION.',
               'REMARK 900 STATIC MODEL ASSEMBLY; NOT AN EXPERIMENTAL COMPLEX.',
               'REMARK 900 UNITS ANGSTROM; RECEPTOR CHAINS A-F; LIGAND CHAIN L, LIG 1.',
               'REMARK 900 RECEPTOR HEAVY ATOMS ONLY; LIGAND RULE-STATE SDF UNCHANGED.',
               'REMARK 900 RECEPTOR SHA256 '+RECEPTOR_PDB_SHA,
               'REMARK 900 LIGAND SDF SHA256 '+sha(source),
               'REMARK 900 EXACT BOND ORDER AND CHARGE REMAIN IN COMPANION SDF.']
        for serial,line in enumerate(receptor,1):lines.append(line[:6]+f'{serial:5d}'+line[11:])
        errors=[];element_counts=Counter()
        for n,atom in enumerate(mol.GetAtoms()):
            pos=conf.GetAtomPosition(n);coords=(pos.x,pos.y,pos.z);serial=15071+n;element=atom.GetSymbol();element_counts[element]+=1
            atomname=(element+str(element_counts[element]))[:4];field=atomname.rjust(4)
            line=f'HETATM{serial:5d} {field} LIG L   1    {coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2}  '
            assert line[21]=='L' and line[17:20]=='LIG'
            published=xyz(line);error=max(abs(x-y) for x,y in zip(coords,published));assert error<=0.0005000001;errors.append(error);lines.append(line)
            atommap.append({'candidate_id':cid,'sdf_atom_1based':n+1,'pdb_serial':serial,'pdb_atom_name':atomname,'element':element,'formal_charge':atom.GetFormalCharge(),'chain':'L','residue_name':'LIG','residue_number':1,'sdf_x_A':coords[0],'sdf_y_A':coords[1],'sdf_z_A':coords[2],'pdb_x_A':published[0],'pdb_y_A':published[1],'pdb_z_A':published[2],'max_coordinate_rounding_error_A':error})
        connectivity={15071+n:[] for n in range(mol.GetNumAtoms())}
        for bond in mol.GetBonds():
            begin,end=bond.GetBeginAtomIdx()+15071,bond.GetEndAtomIdx()+15071
            connectivity[begin].append(end);connectivity[end].append(begin)
            bonds.append({'candidate_id':cid,'sdf_atom1_1based':bond.GetBeginAtomIdx()+1,'sdf_atom2_1based':bond.GetEndAtomIdx()+1,'pdb_serial1':begin,'pdb_serial2':end,'bond_type':str(bond.GetBondType()),'bond_order':bond.GetBondTypeAsDouble(),'aromatic':bond.GetIsAromatic(),'authoritative_bond_source':'ligands/'+source.name})
        for start,neighbours in connectivity.items():
            for offset in range(0,len(neighbours),4):lines.append('CONECT'+f'{start:5d}'+''.join(f'{n:5d}' for n in neighbours[offset:offset+4]))
        lines+=['END']
        output=a.output/'complexes'/(cid+'__M_PUB_ALL__e8__s20260927_complex.pdb')
        output.write_text('\n'.join(lines)+'\n',encoding='ascii')
        index.append({'candidate_id':cid,'source_id':candidate['source_id'],'complex_file':output.relative_to(a.output).as_posix(),'complex_sha256':sha(output),'ligand_sdf_file':sourcecopy.relative_to(a.output).as_posix(),'ligand_sdf_sha256':sha(source),'source_sdf_records':len(supplier),'selected_sdf_record_1based':1,'ligand_source_title':title,'receptor_id':'M_PUB_ALL','receptor_heavy_pdb_sha256':RECEPTOR_PDB_SHA,'receptor_pdbqt_sha256':RECEPTOR_PDBQT_SHA,'receptor_heavy_atoms':15070,'receptor_chains':'A;B;C;D;E;F','ligand_chain':'L','ligand_residue':'LIG 1','ligand_all_atoms':mol.GetNumAtoms(),'ligand_heavy_atoms':mol.GetNumHeavyAtoms(),'ligand_formal_charge':Chem.GetFormalCharge(mol),'coordinate_units':'Angstrom','max_ligand_coordinate_rounding_A':max(errors),'receptor_coordinate_change_A':0,'new_docking':False,'new_prediction':False,'minimization':False,'initial_vina_kcal_mol':candidate['initial_vina_kcal_mol'],'pose_scope':'original first pose, frozen M_PUB_ALL e8 seed20260927','structure_scope':'assembled receptor-plus-ligand static computational model; not experimental structure'})
    writecsv(a.output/'COMPLEX_INDEX.csv',index);writecsv(a.output/'LIGAND_ATOM_MAPPING.csv',atommap);writecsv(a.output/'LIGAND_BOND_MAPPING.csv',bonds)
    receipt={'status':'ASSEMBLED','created_utc':datetime.now(timezone.utc).isoformat(),'complex_count':6,'receptor_source_sha256':RECEPTOR_PDB_SHA,'receptor_pdbqt_source_sha256':RECEPTOR_PDBQT_SHA,'receptor_coordinate_frame_match':True,'receptor_heavy_atoms':15070,'candidate_csv_sha256':sha(a.candidates),'ligand_source_sha256':{r['candidate_id']:r['ligand_sdf_sha256'] for r in index},'script_sha256':sha(Path(__file__)),'software':{'RDKit':rdBase.rdkitVersion},'receptor_coordinate_change_A':0,'max_ligand_rounding_coordinate_A':max(r['max_ligand_coordinate_rounding_A'] for r in index),'new_docking_jobs':0,'new_prediction':False,'minimization':False,'explicit_ligand_hydrogens':'all source SDF atoms retained','bond_orders':'authoritative original SDF copied unchanged; bond mapping also supplied','redistribution_scope':'team-private competition materials; author receptor redistribution license not independently established; not for public GitHub distribution'}
    (a.output/'ASSEMBLY_PROVENANCE.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt,ensure_ascii=False))
if __name__=='__main__':main()
