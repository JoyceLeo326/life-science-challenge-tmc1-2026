
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from chemistry_rules import _heavy_element_bond_query, unspecified_potential_stereo

def choose_new_acyclic_n_geometry(source, state, geometry_index):
    unknown = unspecified_potential_stereo(state)
    match = state.GetSubstructMatch(_heavy_element_bond_query(source), useChirality=False)
    if len(match) != source.GetNumHeavyAtoms():
        raise ValueError('Sensitivity requires unchanged heavy-element/bond graph')
    inverse = {r:s for s,r in enumerate(match)}
    allowed=[]
    for kind,idx in unknown:
        atom=state.GetAtomWithIdx(idx);prior=source.GetAtomWithIdx(inverse[idx])
        carbonyl_adjacent=any(n.GetAtomicNum()==6 and any(b.GetBondType()==Chem.BondType.DOUBLE and b.GetOtherAtom(n).GetAtomicNum() in (8,16) for b in n.GetBonds()) for n in atom.GetNeighbors())
        condition=(kind=='Atom_Tetrahedral' and atom.GetAtomicNum()==7 and atom.GetFormalCharge()==1 and atom.GetDegree()==3 and atom.GetTotalNumHs()==1 and not atom.IsInRing() and not atom.GetIsAromatic() and not carbonyl_adjacent and all(b.GetBondType()==Chem.BondType.SINGLE for b in atom.GetBonds()) and prior.GetAtomicNum()==7 and prior.GetFormalCharge()==0 and prior.GetDegree()==3 and prior.GetTotalNumHs()==0 and prior.GetChiralTag()==Chem.ChiralType.CHI_UNSPECIFIED and not prior.IsInRing())
        if not condition:raise ValueError('Sensitivity cannot relax nonqualifying stable or ring stereochemistry')
        allowed.append(idx)
    if len(allowed)!=1:raise ValueError('Sensitivity preregistered exactly one new acyclic N center')
    options=StereoEnumerationOptions(onlyUnassigned=True,unique=True,tryEmbedding=False,maxIsomers=0)
    variants=sorted(Chem.MolToSmiles(m,isomericSmiles=True) for m in EnumerateStereoisomers(state,options=options))
    if len(variants)!=2:raise ValueError('Expected both and exactly two new N geometry variants')
    chosen=Chem.MolFromSmiles(variants[geometry_index])
    # Clear only newly specified N geometry; all other state stereochemistry must survive.
    query=state.GetSubstructMatch(_heavy_element_bond_query(source),useChirality=False)
    mapping=chosen.GetSubstructMatch(_heavy_element_bond_query(state),useChirality=False)
    cleared=Chem.Mol(chosen)
    for idx in allowed:cleared.GetAtomWithIdx(mapping[idx]).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    if Chem.MolToSmiles(cleared,isomericSmiles=True)!=Chem.MolToSmiles(state,isomericSmiles=True):
        raise ValueError('Sensitivity changed pre-existing stereochemistry or molecular graph')
    return chosen,{'new_N_indices_in_original_rule_state':allowed,'geometry_count':2,'geometry_index':geometry_index,'all_geometry_smiles':variants,'interpretation':'two fixed geometries; no population weighting or inversion-rate claim; excluded from frozen-protocol success rates'}
