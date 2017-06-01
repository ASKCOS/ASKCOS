import numpy as np
import rdkit.Chem as Chem


def score_smiles(smiles):
    x = Chem.MolFromSmiles(smiles)
    total_atoms = x.GetNumHeavyAtoms()
    ring_bonds = sum([b.IsInRing() - b.GetIsAromatic()
                      for b in x.GetBonds()])
    chiral_centers = len(Chem.FindMolChiralCenters(x))

    return - 2.00 * np.power(total_atoms, 1.5) - 1.00 * np.power(ring_bonds, 1.5) - 2.00 * np.power(chiral_centers, 2.0)
