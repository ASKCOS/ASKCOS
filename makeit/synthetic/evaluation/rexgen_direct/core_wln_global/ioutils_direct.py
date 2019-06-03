import rdkit.Chem as Chem
from makeit.synthetic.evaluation.rexgen_direct.core_wln_global.mol_graph import bond_fdim, bond_features
import numpy as np

BOND_TYPE = ["NOBOND", Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE, Chem.rdchem.BondType.AROMATIC] 
N_BOND_CLASS = len(BOND_TYPE)
binary_fdim = 4 + bond_fdim
INVALID_BOND = -1

def get_bin_feature(r, max_natoms):
    comp = {}
    for i, s in enumerate(r.split('.')):
        mol = Chem.MolFromSmiles(s)
        for atom in mol.GetAtoms():
            comp[atom.GetIntProp('molAtomMapNumber') - 1] = i
    n_comp = len(r.split('.'))
    rmol = Chem.MolFromSmiles(r)
    n_atoms = rmol.GetNumAtoms()
    bond_map = {}
    for bond in rmol.GetBonds():
        a1 = bond.GetBeginAtom().GetIntProp('molAtomMapNumber') - 1
        a2 = bond.GetEndAtom().GetIntProp('molAtomMapNumber') - 1
        bond_map[(a1,a2)] = bond_map[(a2,a1)] = bond
        
    features = []
    for i in range(max_natoms):
        for j in range(max_natoms):
            f = np.zeros((binary_fdim,))
            if i >= n_atoms or j >= n_atoms or i == j:
                features.append(f)
                continue
            if (i,j) in bond_map:
                bond = bond_map[(i,j)]
                f[1:1+bond_fdim] = bond_features(bond)
            else:
                f[0] = 1.0
            f[-4] = 1.0 if comp[i] != comp[j] else 0.0
            f[-3] = 1.0 if comp[i] == comp[j] else 0.0
            f[-2] = 1.0 if n_comp == 1 else 0.0
            f[-1] = 1.0 if n_comp > 1 else 0.0
            features.append(f)
    return np.vstack(features).reshape((max_natoms,max_natoms,binary_fdim))

bo_to_index  = {0.0: 0, 1:1, 2:2, 3:3, 1.5:4}
nbos = len(bo_to_index)
def get_bond_label(r, edits, max_natoms):
    rmol = Chem.MolFromSmiles(r)
    n_atoms = rmol.GetNumAtoms()
    rmap = np.zeros((max_natoms, max_natoms, nbos))
    
    for s in edits.split(';'):
        a1,a2,bo = s.split('-')
        x = min(int(a1)-1,int(a2)-1)
        y = max(int(a1)-1, int(a2)-1)
        z = bo_to_index[float(bo)]
        rmap[x,y,z] = rmap[y,x,z] = 1

    labels = []
    sp_labels = []
    for i in range(max_natoms):
        for j in range(max_natoms):
            for k in range(len(bo_to_index)):
                if i == j or i >= n_atoms or j >= n_atoms:
                    labels.append(INVALID_BOND) # mask
                else:
                    labels.append(rmap[i,j,k])
                    if rmap[i,j,k] == 1:
                        sp_labels.append(i * max_natoms * nbos + j * nbos + k)
                        # TODO: check if this is consistent with how TF does flattening
    return np.array(labels), sp_labels

def get_all_batch(re_list):
    mol_list = []
    max_natoms = 0
    for r,e in re_list:
        rmol = Chem.MolFromSmiles(r)
        mol_list.append((r,e))
        if rmol.GetNumAtoms() > max_natoms:
            max_natoms = rmol.GetNumAtoms()
    labels = []
    features = []
    sp_labels = []
    for r,e in mol_list:
        l, sl = get_bond_label(r,e,max_natoms)
        features.append(get_bin_feature(r,max_natoms))
        labels.append(l)
        sp_labels.append(sl)
    return np.array(features), np.array(labels), sp_labels

def get_feature_batch(r_list):
    max_natoms = 0
    for r in r_list:
        rmol = Chem.MolFromSmiles(r)
        if rmol.GetNumAtoms() > max_natoms:
            max_natoms = rmol.GetNumAtoms()

    features = []
    for r in r_list:
        features.append(get_bin_feature(r,max_natoms))
    return np.array(features)
