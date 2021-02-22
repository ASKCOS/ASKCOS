import rdkit.Chem as Chem
from makeit.synthetic.selectivity.mol_graph import bond_fdim, bond_features
import numpy as np

binary_fdim = 5 + bond_fdim

def get_bin_feature(r, max_natoms):
    comp = {}
    rmol = Chem.MolFromSmiles(r)
    for i, frag_ids in enumerate(Chem.GetMolFrags(rmol)):
        for idx in frag_ids:
            comp[idx] = i
    n_comp = len(r.split('.'))
    n_atoms = rmol.GetNumAtoms()
    bond_map = {}
    for bond in rmol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
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

def binary_features_batch(r_list):
    max_natoms = 0
    for r in r_list:
        rmol = Chem.MolFromSmiles(r)
        if rmol.GetNumAtoms() > max_natoms:
            max_natoms = rmol.GetNumAtoms()
    features = []
    for r in r_list:
        features.append(get_bin_feature(r,max_natoms))
    return np.array(features)
