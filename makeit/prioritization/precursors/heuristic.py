from makeit.prioritization.prioritizer import Prioritizer
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from makeit.utilities.io.logging import MyLogger
heuristic_precursor_prioritizer_loc = 'heuristic_precursor_prioritizer'


class HeuristicPrecursorPrioritizer(Prioritizer):

    def __init__(self):
        MyLogger.print_and_log(
            'Using heuristic prioritization method for directing the tree expansion.', heuristic_precursor_prioritizer_loc)

    def get_priority(self, retroPrecursor):

        necessary_reagent_atoms = retroPrecursor.necessary_reagent.count(
            '[') / 2.
        scores = []
        for smiles in retroPrecursor.smiles_list:
            x = Chem.MolFromSmiles(smiles)
            total_atoms = x.GetNumHeavyAtoms()
            ring_bonds = sum([b.IsInRing() - b.GetIsAromatic()
                              for b in x.GetBonds()])
            chiral_centers = len(Chem.FindMolChiralCenters(x))

            scores.append(
                - 2.00 * np.power(total_atoms, 1.5)
                - 1.00 * np.power(ring_bonds, 1.5)
                - 2.00 * np.power(chiral_centers, 2.0)
            )

        return np.min(scores) - 4.00 * np.power(necessary_reagent_atoms, 2.0)

    def load_model(self):
        pass
