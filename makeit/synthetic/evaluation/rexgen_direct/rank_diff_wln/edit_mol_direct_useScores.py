import rdkit
from rdkit import Chem
from optparse import OptionParser

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(4)

BOND_TYPE = [0, Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE, Chem.rdchem.BondType.AROMATIC] 
BOND_FLOAT_TO_TYPE = {
    0.0: BOND_TYPE[0],
    1.0: BOND_TYPE[1],
    2.0: BOND_TYPE[2],
    3.0: BOND_TYPE[3],
    1.5: BOND_TYPE[4],
}
def copy_edit_mol(mol):
    new_mol = Chem.RWMol(Chem.MolFromSmiles(''))
    for atom in mol.GetAtoms():
        new_atom = Chem.Atom(atom.GetSymbol())
        new_atom.SetFormalCharge(atom.GetFormalCharge()) # TODO: How to deal with changing formal charge?
        new_atom.SetAtomMapNum(atom.GetAtomMapNum())
        new_mol.AddAtom(new_atom)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        bt = bond.GetBondType()
        new_mol.AddBond(a1, a2, bt)
    return new_mol

def get_product_smiles(rmol, edits, tatoms):
    smiles = edit_mol(rmol, edits, tatoms)
    if len(smiles) != 0: return smiles
    try:
        Chem.Kekulize(rmol)
    except Exception as e:
        return smiles
    return edit_mol(rmol, edits, tatoms)

def edit_mol(rmol, edits, tatoms):
    #new_mol = copy_edit_mol(rmol)
    new_mol = Chem.RWMol(rmol)
    [a.SetNumExplicitHs(0) for a in new_mol.GetAtoms()]

    amap = {}
    for atom in rmol.GetAtoms():
        amap[atom.GetAtomMapNum() - 1] = atom.GetIdx()

    for x,y,t,v in edits:
        bond = new_mol.GetBondBetweenAtoms(amap[x],amap[y])
        # a1 = new_mol.GetAtomWithIdx(amap[x])
        # a2 = new_mol.GetAtomWithIdx(amap[y])
        if bond is not None:
            new_mol.RemoveBond(amap[x],amap[y])
        if t > 0:
            new_mol.AddBond(amap[x],amap[y],BOND_FLOAT_TO_TYPE[t])
    
    pred_mol = new_mol.GetMol()
    pred_smiles = Chem.MolToSmiles(pred_mol)
    pred_list = pred_smiles.split('.')
    pred_mols = []
    for pred_smiles in pred_list:
        mol = Chem.MolFromSmiles(pred_smiles)
        if mol is None: continue
        atom_set = set([atom.GetAtomMapNum() - 1 for atom in mol.GetAtoms()])
        if len(atom_set & tatoms) == 0:
            continue
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        pred_mols.append(mol)

    return '.'.join( sorted([Chem.MolToSmiles(pred_mol) for pred_mol in pred_mols]) )
