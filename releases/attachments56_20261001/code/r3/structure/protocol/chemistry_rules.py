"""Narrow rule for transient N+H geometry after protonating a tertiary amine.

This is a computational single-geometry approximation, not a microstate model.
The rule does not relax source-defined carbon or stable quaternary/bridgehead N.
"""
from rdkit import Chem

RULE_VERSION = "exchangeable_tertiary_ring_NH_and_full_potential_stereo_v3_20260929"


def unspecified_potential_stereo(mol):
    """Report atom and double-bond stereogenic elements lacking source assignment."""
    return [(str(info.type), int(info.centeredOn)) for info in Chem.FindPotentialStereo(mol)
            if str(info.specified) == "Unspecified"]


def _heavy_element_bond_query(mol):
    query = Chem.RWMol()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        query.AddAtom(Chem.AtomFromSmarts(f"[#{atom.GetAtomicNum()}]"))
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if mol.GetAtomWithIdx(i).GetAtomicNum() == 1 or mol.GetAtomWithIdx(j).GetAtomicNum() == 1:
            continue
        query.AddBond(i, j, bond.GetBondType())
    return query.GetMol()


def exchangeable_n_indices(source, rule):
    """Return rule-state atom indices meeting the verified narrow N+H condition.

    Full heavy-element/bond-order graph must map to the neutral input; protonation
    is the only permitted change at the candidate N. Ambiguous source stereo is
    checked separately by callers.
    """
    if source.GetNumHeavyAtoms() != rule.GetNumHeavyAtoms():
        return []
    match = rule.GetSubstructMatch(_heavy_element_bond_query(source), useChirality=False)
    if len(match) != source.GetNumHeavyAtoms():
        return []
    source_for_rule = {rule_idx: source_idx for source_idx, rule_idx in enumerate(match)}
    rings = [set(r) for r in Chem.GetSymmSSSR(rule)]
    allowed = []
    for atom in rule.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1:
            continue
        if atom.GetDegree() != 3 or atom.GetTotalNumHs() != 1:
            continue
        if sum(idx in ring for ring in rings) != 1:
            continue
        if any(b.GetBondType() != Chem.BondType.SINGLE for b in atom.GetBonds()):
            continue
        prior = source.GetAtomWithIdx(source_for_rule[idx])
        if prior.GetAtomicNum() != 7 or prior.GetFormalCharge() != 0:
            continue
        if prior.GetDegree() != 3 or prior.GetTotalNumHs() != 0:
            continue
        if any(b.GetBondType() != Chem.BondType.SINGLE for b in prior.GetBonds()):
            continue
        allowed.append(idx)
    return allowed


def stable_identity_smiles(mol, exchangeable_indices):
    copy = Chem.Mol(mol)
    for idx in exchangeable_indices:
        copy.GetAtomWithIdx(idx).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    Chem.AssignStereochemistry(copy, cleanIt=True, force=True)
    return Chem.MolToSmiles(copy, isomericSmiles=True)
