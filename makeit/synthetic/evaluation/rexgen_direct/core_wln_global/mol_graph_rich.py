import rdkit
import rdkit.Chem as Chem
import numpy as np

# Extra atom features
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
import rdkit.Chem.EState as EState
import rdkit.Chem.rdPartialCharges as rdPartialCharges
import rdkit.Chem.rdChemReactions as rdRxns

elem_list = ['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na', 'Ca', 'Fe', 'As', 'Al', 'I', 'B', 'V', 'K', 'Tl', 'Yb', 'Sb', 'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti', 'Zn', 'H', 'Li', 'Ge', 'Cu', 'Au', 'Ni', 'Cd', 'In', 'Mn', 'Zr', 'Cr', 'Pt', 'Hg', 'Pb', 'W', 'Ru', 'Nb', 'Re', 'Te', 'Rh', 'Tc', 'Ba', 'Bi', 'Hf', 'Mo', 'U', 'Sm', 'Os', 'Ir', 'Ce','Gd','Ga','Cs', 'unknown']
max_nb = 10

def onek_encoding_unk(x, allowable_set):
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))

def atom_features(atom):
    attributes = onek_encoding_unk(atom.GetSymbol(), elem_list) \
            + onek_encoding_unk(atom.GetDegree(), [0,1,2,3,4,5]) \
            + onek_encoding_unk(atom.GetExplicitValence(), [1,2,3,4,5,6]) \
            + onek_encoding_unk(atom.GetImplicitValence(), [0,1,2,3,4,5]) \
            + [atom.GetIsAromatic()] \
            + [atom.GetIsAromatic() == False and any([neighbor.GetIsAromatic() for neighbor in atom.GetNeighbors()])] \
            + [atom.IsInRing()] \
            + [atom.GetAtomicNum() in [9, 17, 35, 53, 85, 117]] \
            + [atom.GetAtomicNum() in [8, 16, 34, 52, 84, 116]] \
            + [atom.GetAtomicNum() in [7, 15, 33, 51, 83]] \
            + [atom.GetAtomicNum() in [3, 11, 19, 37, 55, 87]] \
            + [atom.GetAtomicNum() in [4, 12, 20, 38, 56, 88]] \
            + [atom.GetAtomicNum() in [13, 22, 24, 25, 26, 27, 28, 29, 30, 33, 42, 44, 45, 46, 47, 48, 49, 50, 78, 80, 82]] \
            + [atom.GetDoubleProp('crippen_logp'), atom.GetDoubleProp('crippen_mr'), atom.GetDoubleProp('tpsa'),
               atom.GetDoubleProp('asa'), atom.GetDoubleProp('estate'),
               atom.GetDoubleProp('_GasteigerCharge'), atom.GetDoubleProp('_GasteigerHCharge')]
    attributes = np.array(attributes, dtype=np.float32)
    attributes[np.isnan(attributes)] = 0.0 # filter nan
    attributes[np.isinf(attributes)] = 0.0 # filter inf
    return attributes

def assignProperties(mol):
    '''
    Calculate atom-level descriptors that can be used in featurization
    '''
    for (i, x) in enumerate(rdMolDescriptors._CalcCrippenContribs(mol)):
        mol.GetAtomWithIdx(i).SetDoubleProp('crippen_logp',x[0])
        mol.GetAtomWithIdx(i).SetDoubleProp('crippen_mr', x[1])
    for (i, x) in enumerate(rdMolDescriptors._CalcTPSAContribs(mol)):
        mol.GetAtomWithIdx(i).SetDoubleProp('tpsa', x)
    for (i, x) in enumerate(rdMolDescriptors._CalcLabuteASAContribs(mol)[0]):
        mol.GetAtomWithIdx(i).SetDoubleProp('asa', x)
    for (i, x) in enumerate(EState.EStateIndices(mol)):
        mol.GetAtomWithIdx(i).SetDoubleProp('estate', x)
    rdPartialCharges.ComputeGasteigerCharges(mol) # '_GasteigerCharge', '_GasteigerHCharge'



def bond_features(bond):
    bt = bond.GetBondType()
    return np.array([bt == Chem.rdchem.BondType.SINGLE, bt == Chem.rdchem.BondType.DOUBLE,
                     bt == Chem.rdchem.BondType.TRIPLE, bt == Chem.rdchem.BondType.AROMATIC,
                     bond.GetIsConjugated(), bond.IsInRing()], dtype=np.float32)

def smiles2graph(smiles, idxfunc=lambda x:x.GetIdx()):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("Could not parse smiles string:", smiles)

    n_atoms = mol.GetNumAtoms()
    n_bonds = max(mol.GetNumBonds(), 1)
    fatoms = np.zeros((n_atoms, atom_fdim))
    fbonds = np.zeros((n_bonds, bond_fdim))
    atom_nb = np.zeros((n_atoms, max_nb), dtype=np.int32)
    bond_nb = np.zeros((n_atoms, max_nb), dtype=np.int32)
    num_nbs = np.zeros((n_atoms,), dtype=np.int32)

    assignProperties(mol)

    for atom in mol.GetAtoms():
        idx = idxfunc(atom)
        if idx >= n_atoms:
            raise Exception(smiles)
        fatoms[idx] = atom_features(atom)

    for bond in mol.GetBonds():
        a1 = idxfunc(bond.GetBeginAtom())
        a2 = idxfunc(bond.GetEndAtom())
        idx = bond.GetIdx()
        if num_nbs[a1] == max_nb or num_nbs[a2] == max_nb:
            raise Exception(smiles)
        atom_nb[a1,num_nbs[a1]] = a2
        atom_nb[a2,num_nbs[a2]] = a1
        bond_nb[a1,num_nbs[a1]] = idx
        bond_nb[a2,num_nbs[a2]] = idx
        num_nbs[a1] += 1
        num_nbs[a2] += 1
        fbonds[idx] = bond_features(bond)
    return fatoms, fbonds, atom_nb, bond_nb, num_nbs

def pack2D(arr_list):
    N = max([x.shape[0] for x in arr_list])
    M = max([x.shape[1] for x in arr_list])
    a = np.zeros((len(arr_list), N, M))
    for i, arr in enumerate(arr_list):
        n = arr.shape[0]
        m = arr.shape[1]
        a[i,0:n,0:m] = arr
    return a

def pack2D_withidx(arr_list):
    N = max([x.shape[0] for x in arr_list])
    M = max([x.shape[1] for x in arr_list])
    a = np.zeros((len(arr_list), N, M, 2))
    for i, arr in enumerate(arr_list):
        n = arr.shape[0]
        m = arr.shape[1]
        a[i,0:n,0:m,0] = i
        a[i,0:n,0:m,1] = arr
    return a

def pack1D(arr_list):
    N = max([x.shape[0] for x in arr_list])
    a = np.zeros((len(arr_list), N))
    for i, arr in enumerate(arr_list):
        n = arr.shape[0]
        a[i,0:n] = arr
    return a

def get_mask(arr_list):
    N = max([x.shape[0] for x in arr_list])
    a = np.zeros((len(arr_list), N))
    for i, arr in enumerate(arr_list):
        for j in range(arr.shape[0]):
            a[i][j] = 1
    return a

def smiles2graph_list(smiles_list, idxfunc=lambda x:x.GetIdx()):
    res = list(map(lambda x:smiles2graph(x,idxfunc), smiles_list))
    fatom_list, fbond_list, gatom_list, gbond_list, nb_list = zip(*res)
    return pack2D(fatom_list), pack2D(fbond_list), pack2D_withidx(gatom_list), pack2D_withidx(gbond_list), pack1D(nb_list), get_mask(fatom_list)

m = Chem.MolFromSmiles('CC')
assignProperties(m)
atom = m.GetAtoms()[0]
bond = m.GetBonds()[0]
atom_fdim = len(atom_features(atom))
bond_fdim = len(bond_features(bond))

if __name__ == "__main__":
    np.set_printoptions(threshold='nan')
    a,b,c,d,e,f = smiles2graph_list(["c1cccnc1",'c1nccc2n1ccc2'])
    print(a)
    print(b)
    print(c)
    print(d)
    print(e)
    print(f)
