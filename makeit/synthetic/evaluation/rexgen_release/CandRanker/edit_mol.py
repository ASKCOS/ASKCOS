import rdkit
from rdkit import Chem
from optparse import OptionParser

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(4)

BOND_TYPE = [0, Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE, Chem.rdchem.BondType.AROMATIC] 

def edit_mol(rmol, edits):
    n_atoms = rmol.GetNumAtoms()

    new_mol = Chem.RWMol(rmol)
    amap = {}
    numH = {}
    for atom in rmol.GetAtoms():
        amap[atom.GetIntProp('molAtomMapNumber')] = atom.GetIdx()
        numH[atom.GetIntProp('molAtomMapNumber')] = atom.GetNumExplicitHs()

    for x,y,t in edits:
        bond = new_mol.GetBondBetweenAtoms(amap[x],amap[y])
        a1 = new_mol.GetAtomWithIdx(amap[x])
        a2 = new_mol.GetAtomWithIdx(amap[y])
        if bond is not None:
            val = BOND_TYPE.index(bond.GetBondType())
            new_mol.RemoveBond(amap[x],amap[y])
            numH[x] += val
            numH[y] += val

        if t > 0:
            new_mol.AddBond(amap[x],amap[y],BOND_TYPE[t])
            numH[x] -= t
            numH[y] -= t

    for atom in new_mol.GetAtoms():
        val = numH[atom.GetIntProp('molAtomMapNumber')]
        if val >= 0:
            atom.SetNumExplicitHs(val)

    pred_mol = new_mol.GetMol()
    for atom in pred_mol.GetAtoms():
        atom.ClearProp('molAtomMapNumber')
    pred_smiles = Chem.MolToSmiles(pred_mol)

    return pred_smiles

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-t", "--pred", dest="pred_path")
    parser.add_option("-g", "--gold", dest="gold_path")
    opts,args = parser.parse_args()

    fpred = open(opts.pred_path)
    fgold = open(opts.gold_path)

    rank = []
    for line in fpred:
        line = line.strip('\r\n |')
        gold = fgold.readline()
        rex,_ = gold.split()
        r,_,p = rex.split('>')
        rmol = Chem.MolFromSmiles(r)
        pmol = Chem.MolFromSmiles(p)
        
        patoms = set()
        pbonds = {}
        for bond in pmol.GetBonds():
            a1 = bond.GetBeginAtom().GetIntProp('molAtomMapNumber')
            a2 = bond.GetEndAtom().GetIntProp('molAtomMapNumber') 
            t = BOND_TYPE.index(bond.GetBondType())
            a1,a2 = min(a1,a2),max(a1,a2)
            pbonds[(a1,a2)] = t
            patoms.add(a1)
            patoms.add(a2)

        rbonds = {}
        for bond in rmol.GetBonds():
            a1 = bond.GetBeginAtom().GetIntProp('molAtomMapNumber')
            a2 = bond.GetEndAtom().GetIntProp('molAtomMapNumber') 
            t = BOND_TYPE.index(bond.GetBondType())
            a1,a2 = min(a1,a2),max(a1,a2)
            if a1 in patoms or a2 in patoms:
                rbonds[(a1,a2)] = t
        
        rk = 10
        for idx,edits in enumerate(line.split('|')):
            cbonds = [] 
            pred = dict(rbonds)
            for edit in edits.split():
                x,y,t = edit.split('-')
                x,y,t = int(x),int(y),int(t)
                cbonds.append((x,y,t))
                if t == 0 and (x,y) in rbonds:
                    del pred[(x,y)]
                if t > 0:
                    pred[(x,y)] = t

            if pred == pbonds: 
                rk = idx + 1
                break
            for atom in pmol.GetAtoms():
                atom.ClearProp('molAtomMapNumber')
            psmiles = Chem.MolToSmiles(pmol)
            psmiles = set(psmiles.split('.'))
            pred_smiles = set(edit_mol(r, cbonds).split('.'))
            if psmiles <= pred_smiles:
                rk = idx + 1
                break
        rank.append(rk)

    n = 1.0 * len(rank)
    top1,top3,top5 = 0,0,0
    for idx in rank:
        if idx == 1: top1 += 1
        if idx <= 3: top3 += 1
        if idx <= 5: top5 += 1

    print('%.4f, %.4f, %.4f' % (top1 / n, top3 / n, top5 / n))
