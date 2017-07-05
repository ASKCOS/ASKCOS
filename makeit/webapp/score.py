import numpy as np
import rdkit.Chem as Chem


def score_smiles(smiles, Pricer=None):
    if Pricer is not None:
        ppg = Pricer.lookup_smiles(smiles, alreadyCanonical=True)
        if ppg:
            return - ppg / 5.0 # basically "free"
    x = Chem.MolFromSmiles(smiles)
    total_atoms = x.GetNumHeavyAtoms()
    ring_bonds = sum([b.IsInRing() - b.GetIsAromatic()
                      for b in x.GetBonds()])
    chiral_centers = len(Chem.FindMolChiralCenters(x))

    return - 2.00 * np.power(total_atoms, 1.5) - 1.00 * np.power(ring_bonds, 1.5) - 2.00 * np.power(chiral_centers, 2.0)
