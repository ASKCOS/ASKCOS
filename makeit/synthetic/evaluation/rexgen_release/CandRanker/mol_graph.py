import rdkit
import rdkit.Chem as Chem
import numpy as np
import random

elem_list = ['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na', 'Ca', 'Fe', 'As', 'Al', 'I', 'B', 'V', 'K', 'Tl', 'Yb', 'Sb', 'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti', 'Zn', 'H', 'Li', 'Ge', 'Cu', 'Au', 'Ni', 'Cd', 'In', 'Mn', 'Zr', 'Cr', 'Pt', 'Hg', 'Pb', 'W', 'Ru', 'Nb', 'Re', 'Te', 'Rh', 'Tc', 'Ba', 'Bi', 'Hf', 'Mo', 'U', 'Sm', 'Os', 'Ir', 'Ce','Gd','Ga','Cs', 'unknown']
bond_types = [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE, Chem.rdchem.BondType.AROMATIC]
atom_fdim = len(elem_list) + 6 + 6 + 6 + 1
bond_fdim = 5
max_nb = 10

def onek_encoding_unk(x, allowable_set):
    if x not in allowable_set:
        x = allowable_set[-1]
    return map(lambda s: x == s, allowable_set)

def atom_features(atom):
    return np.array(onek_encoding_unk(atom.GetSymbol(), elem_list) 
            + onek_encoding_unk(atom.GetDegree(), [0,1,2,3,4,5]) 
            + onek_encoding_unk(atom.GetExplicitValence(), [1,2,3,4,5,6])
            + onek_encoding_unk(atom.GetImplicitValence(), [0,1,2,3,4,5])
            + [atom.GetIsAromatic()], dtype=np.float32)

def bond_features(bond):
    bt = bond.GetBondType()
    return np.array([bt == Chem.rdchem.BondType.SINGLE, bt == Chem.rdchem.BondType.DOUBLE, bt == Chem.rdchem.BondType.TRIPLE, bt == Chem.rdchem.BondType.AROMATIC, bond.IsInRing()], dtype=np.float32)

def search(buf, cur_bonds, core_bonds, free_vals, depth):
    if depth >= len(core_bonds):
        buf.append([u for u in cur_bonds])
        return
    x,y,t = core_bonds[depth]
    if t >= 0:
        cur_bonds.append((x,y,t))
        free_vals[x] -= t
        free_vals[y] -= t
        search(buf, cur_bonds, core_bonds, free_vals, depth + 1)
        free_vals[x] += t
        free_vals[y] += t
        cur_bonds.pop()
    else:
        for k in xrange(4):
            if k > free_vals[x] or k > free_vals[y]:
                break
            cur_bonds.append((x,y,k))
            free_vals[x] -= k
            free_vals[y] -= k
            search(buf, cur_bonds, core_bonds, free_vals, depth + 1)
            free_vals[x] += k
            free_vals[y] += k
            cur_bonds.pop()
            
def packnb(arr_list):
    N = max([x.shape[0] for x in arr_list])
    M = max([x.shape[1] for x in arr_list])
    a = np.zeros((len(arr_list), N, M, 2))
    for i, arr in enumerate(arr_list):
        n = arr.shape[0]
        m = arr.shape[1]
        a[i,0:n,0:m,0] = i
        a[i,0:n,0:m,1] = arr
    return a

def floodfill(cur_id, cur_label, comp, core_bonds):
    comp[cur_id] = cur_label
    x,y = core_bonds[cur_id]
    for i in xrange(len(core_bonds)):
        if comp[i] >= 0: continue
        u,v = core_bonds[i]
        if x == u or x == v or y == u or y == v:
            floodfill(i, cur_label, comp, core_bonds)

def smiles2graph(rsmiles, psmiles, core_bonds, cutoff=500, idxfunc=lambda x:x.GetIntProp('molAtomMapNumber') - 1):
    mol = Chem.MolFromSmiles(rsmiles)
    pmol = Chem.MolFromSmiles(psmiles)
    if not mol or not pmol:
        raise ValueError("Could not parse smiles string:", rsmiles + '>>' + psmiles)

    n_atoms = mol.GetNumAtoms()
    n_bonds = max(mol.GetNumBonds(), 1)
    fatoms = np.zeros((n_atoms, atom_fdim))
    fbonds = np.zeros((n_bonds, bond_fdim))
    atom_nb = np.zeros((n_atoms, max_nb), dtype=np.int32)
    bond_nb = np.zeros((n_atoms, max_nb), dtype=np.int32)
    num_nbs = np.zeros((n_atoms,), dtype=np.int32)
    raw_atom_nb = np.zeros((n_atoms, max_nb), dtype=np.int32)
    raw_bond_nb = np.zeros((n_atoms, max_nb), dtype=np.int32)
    raw_num_nbs = np.zeros((n_atoms,), dtype=np.int32)
    free_vals = np.zeros((n_atoms,))
    pfree_vals = np.zeros((n_atoms,))

    gbonds = {(x,y):0 for x,y in core_bonds}

    #Feature Extraction
    for atom in mol.GetAtoms():
        idx = idxfunc(atom)
        fatoms[idx] = atom_features(atom)
        free_vals[idx] += atom.GetTotalNumHs() + abs(atom.GetFormalCharge())
    
    tatoms = set()
    #Calculate free slots for each atom in product
    for bond in pmol.GetBonds():
        a1 = idxfunc(bond.GetBeginAtom())
        a2 = idxfunc(bond.GetEndAtom())
        t = bond_types.index(bond.GetBondType()) + 1
        a1,a2 = min(a1,a2),max(a1,a2)
        tatoms.add(a1)
        tatoms.add(a2)
        if (a1,a2) in core_bonds:
            gbonds[(a1,a2)] = t
            tval = t if t < 4 else 1.5
            pfree_vals[a1] += tval
            pfree_vals[a2] += tval
    
    rbonds = {}
    ring_bonds = set()
    #Calculate free slots for each atom in reactant
    for bond in mol.GetBonds():
        idx = bond.GetIdx()
        a1 = idxfunc(bond.GetBeginAtom())
        a2 = idxfunc(bond.GetEndAtom())
        t = bond_types.index(bond.GetBondType())
        a1,a2 = min(a1,a2),max(a1,a2)
        tval = t + 1 if t < 3 else 1.5
        rbonds[(a1,a2)] = t + 1
        if (a1,a2) in core_bonds:
            free_vals[a1] += tval 
            free_vals[a2] += tval
        if bond.IsInRing():
            ring_bonds.add((a1,a2))

    #Fix golden label
    for x,y in gbonds.iterkeys():
        if x not in tatoms and y not in tatoms and (x,y) in rbonds:
            gbonds[(x,y)] = rbonds[(x,y)]

    #Take the max just in case
    for i in xrange(n_atoms):
        if pfree_vals[i] > free_vals[i]:
            free_vals[i] = pfree_vals[i]

    #Gather the golden label
    gold_bonds = set()
    for x,y in core_bonds:
        t = gbonds[(x,y)]
        gold_bonds.add( (x,y,t) )

    #Get connected components in core bonds
    comp = [-1] * len(core_bonds)
    tot = 0
    for i in xrange(len(core_bonds)):
        if comp[i] == -1:
            floodfill(i, tot, comp, core_bonds)
            tot += 1
    
    core_configs = []
    for cur_id in xrange(tot):
        cand_bonds = []
        for i in xrange(len(core_bonds)):
            x,y = core_bonds[i]
            if comp[i] == cur_id: t = -1
            elif (x,y) not in rbonds: t = 0
            else: t = rbonds[(x,y)]
            cand_bonds.append((x,y,t))
        search(core_configs, [], cand_bonds, free_vals, 0)
    
    random.shuffle(core_configs)
    idx = -1
    for i,cand_bonds in enumerate(core_configs):
        if set(cand_bonds) == gold_bonds:
            idx = i
            break

    #Reaction Core Miss
    if idx == -1:
        core_configs = [list(gold_bonds)] + core_configs
    else:
        core_configs[0], core_configs[idx] = core_configs[idx], core_configs[0]

    if cutoff > 0:
        core_configs = core_configs[:cutoff]

    n_batch = len(core_configs) + 1
    labels = np.zeros((n_batch-1,))
    labels[0] = 1

    for bond in mol.GetBonds():
        idx = bond.GetIdx()
        a1 = idxfunc(bond.GetBeginAtom())
        a2 = idxfunc(bond.GetEndAtom())
        a1,a2 = min(a1,a2),max(a1,a2)

        if (a1,a2) not in core_bonds:
            raw_atom_nb[a1,raw_num_nbs[a1]] = a2
            raw_atom_nb[a2,raw_num_nbs[a2]] = a1
            raw_bond_nb[a1,raw_num_nbs[a1]] = idx
            raw_bond_nb[a2,raw_num_nbs[a2]] = idx
            raw_num_nbs[a1] += 1
            raw_num_nbs[a2] += 1

        atom_nb[a1,num_nbs[a1]] = a2
        atom_nb[a2,num_nbs[a2]] = a1
        bond_nb[a1,num_nbs[a1]] = idx
        bond_nb[a2,num_nbs[a2]] = idx
        num_nbs[a1] += 1
        num_nbs[a2] += 1
        fbonds[idx] = bond_features(bond)
    
    num_newbonds = len(core_bonds) * 2 + 1
    new_fbonds = np.zeros((n_bonds+num_newbonds, bond_fdim))
    new_fbonds[:n_bonds,:] = fbonds
    fbonds = new_fbonds
    batch_fbonds, batch_anb, batch_bnb, batch_nbs = [fbonds], [atom_nb], [bond_nb], [num_nbs]

    for core_bonds in core_configs:
        atom_nb2 = np.copy(raw_atom_nb)
        bond_nb2 = np.copy(raw_bond_nb)
        num_nbs2 = np.copy(raw_num_nbs)
        fbonds2 = np.copy(fbonds)
        n_bonds2 = n_bonds + 1
        
        for x,y,t in core_bonds:
            if t == 0: continue
            atom_nb2[x,num_nbs2[x]] = y
            atom_nb2[y,num_nbs2[y]] = x
            bond_nb2[x,num_nbs2[x]] = n_bonds2
            bond_nb2[y,num_nbs2[y]] = n_bonds2
            num_nbs2[x] += 1
            num_nbs2[y] += 1
            fbonds2[n_bonds2] = onek_encoding_unk(t-1, [0,1,2,3,4]) 
            if (x,y) in ring_bonds:
                fbonds2[n_bonds2][4] = 1
            n_bonds2 += 1
        batch_fbonds.append(fbonds2)
        batch_anb.append(atom_nb2)
        batch_bnb.append(bond_nb2)
        batch_nbs.append(num_nbs2)

    return (np.array([fatoms] * n_batch), np.array(batch_fbonds), packnb(batch_anb), packnb(batch_bnb), np.array(batch_nbs), labels), core_configs

def smiles2graph_test(rsmiles, core_bonds, idxfunc=lambda x:x.GetIntProp('molAtomMapNumber') - 1):
    mol = Chem.MolFromSmiles(rsmiles)

    n_atoms = mol.GetNumAtoms()
    n_bonds = max(mol.GetNumBonds(), 1)
    fatoms = np.zeros((n_atoms, atom_fdim))
    fbonds = np.zeros((n_bonds, bond_fdim))
    atom_nb = np.zeros((n_atoms, max_nb), dtype=np.int32)
    bond_nb = np.zeros((n_atoms, max_nb), dtype=np.int32)
    num_nbs = np.zeros((n_atoms,), dtype=np.int32)
    raw_atom_nb = np.zeros((n_atoms, max_nb), dtype=np.int32)
    raw_bond_nb = np.zeros((n_atoms, max_nb), dtype=np.int32)
    raw_num_nbs = np.zeros((n_atoms,), dtype=np.int32)
    free_vals = np.zeros((n_atoms,))

    #Feature Extraction
    for atom in mol.GetAtoms():
        idx = idxfunc(atom)
        fatoms[idx] = atom_features(atom)
        free_vals[idx] += atom.GetTotalNumHs() + abs(atom.GetFormalCharge())
    
    rbonds = {}
    ring_bonds = set()
    #Calculate free slots for each atom in reactant
    for bond in mol.GetBonds():
        idx = bond.GetIdx()
        a1 = idxfunc(bond.GetBeginAtom())
        a2 = idxfunc(bond.GetEndAtom())
        t = bond_types.index(bond.GetBondType())
        a1,a2 = min(a1,a2),max(a1,a2)
        tval = t + 1 if t < 3 else 1.5
        rbonds[(a1,a2)] = t + 1
        if (a1,a2) in core_bonds:
            free_vals[a1] += tval 
            free_vals[a2] += tval
        if bond.IsInRing():
            ring_bonds.add((a1,a2))

    #Get connected components in core bonds
    comp = [-1] * len(core_bonds)
    tot = 0
    for i in xrange(len(core_bonds)):
        if comp[i] == -1:
            floodfill(i, tot, comp, core_bonds)
            tot += 1
    
    core_configs = []
    for cur_id in xrange(tot):
        cand_bonds = []
        for i in xrange(len(core_bonds)):
            x,y = core_bonds[i]
            if comp[i] == cur_id: t = -1
            elif (x,y) not in rbonds: t = 0
            else: t = rbonds[(x,y)]
            cand_bonds.append((x,y,t))
        search(core_configs, [], cand_bonds, free_vals, 0)

    n_batch = len(core_configs) + 1

    for bond in mol.GetBonds():
        idx = bond.GetIdx()
        a1 = idxfunc(bond.GetBeginAtom())
        a2 = idxfunc(bond.GetEndAtom())
        a1,a2 = min(a1,a2),max(a1,a2)

        if (a1,a2) not in core_bonds:
            raw_atom_nb[a1,raw_num_nbs[a1]] = a2
            raw_atom_nb[a2,raw_num_nbs[a2]] = a1
            raw_bond_nb[a1,raw_num_nbs[a1]] = idx
            raw_bond_nb[a2,raw_num_nbs[a2]] = idx
            raw_num_nbs[a1] += 1
            raw_num_nbs[a2] += 1

        atom_nb[a1,num_nbs[a1]] = a2
        atom_nb[a2,num_nbs[a2]] = a1
        bond_nb[a1,num_nbs[a1]] = idx
        bond_nb[a2,num_nbs[a2]] = idx
        num_nbs[a1] += 1
        num_nbs[a2] += 1
        fbonds[idx] = bond_features(bond)
    
    num_newbonds = len(core_bonds) * 2 + 1
    new_fbonds = np.zeros((n_bonds+num_newbonds, bond_fdim))
    new_fbonds[:n_bonds,:] = fbonds
    fbonds = new_fbonds
    batch_fbonds, batch_anb, batch_bnb, batch_nbs = [fbonds], [atom_nb], [bond_nb], [num_nbs]

    for core_bonds in core_configs:
        atom_nb2 = np.copy(raw_atom_nb)
        bond_nb2 = np.copy(raw_bond_nb)
        num_nbs2 = np.copy(raw_num_nbs)
        fbonds2 = np.copy(fbonds)
        n_bonds2 = n_bonds + 1
        
        for x,y,t in core_bonds:
            if t == 0: continue
            atom_nb2[x,num_nbs2[x]] = y
            atom_nb2[y,num_nbs2[y]] = x
            bond_nb2[x,num_nbs2[x]] = n_bonds2
            bond_nb2[y,num_nbs2[y]] = n_bonds2
            num_nbs2[x] += 1
            num_nbs2[y] += 1
            fbonds2[n_bonds2] = onek_encoding_unk(t-1, [0,1,2,3,4]) 
            if (x,y) in ring_bonds:
                fbonds2[n_bonds2][4] = 1
            n_bonds2 += 1
        batch_fbonds.append(fbonds2)
        batch_anb.append(atom_nb2)
        batch_bnb.append(bond_nb2)
        batch_nbs.append(num_nbs2)

    return (np.array([fatoms] * n_batch), np.array(batch_fbonds), packnb(batch_anb), packnb(batch_bnb), np.array(batch_nbs)), core_configs

if __name__ == "__main__":
    a,b,c,d,e,f = smiles2graph("[OH:1][CH3:2]", "[O:1]=[CH2:2]", [(0,1)])[0]
    print a
    print b
    print c
    print d
    print e
    print f
