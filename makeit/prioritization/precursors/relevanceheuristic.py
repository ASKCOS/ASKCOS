from makeit.prioritization.prioritizer import Prioritizer
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from makeit.utilities.buyable.pricer import Pricer
from makeit.utilities.io.logger import MyLogger
heuristic_precursor_prioritizer_loc = 'relevanceheuristic_precursor_prioritizer'


class RelevanceHeuristicPrecursorPrioritizer(Prioritizer):
    """A precursor Prioritizer that uses a heuristic and template relevance.

    Attributes:
        pricer (Pricer or None): Used to look up chemical prices.
    """
    def __init__(self):
        """Initializes RelevanceHeuristicPrecursorPrioritizer."""
        self.pricer = None
        self._loaded = False

    def score_precursor(self, precursor):
        scores = []
        necessary_reagent_atoms = precursor['necessary_reagent'].count('[')/2.
        for smiles in precursor['smiles_list']:
            ppg = self.pricer.lookup_smiles(smiles, alreadyCanonical=True)
            # If buyable, basically free
            if ppg:
                scores.append(- ppg / 1000.0)
                continue

            # Else, use heuristic
            mol = Chem.MolFromSmiles(smiles)
            total_atoms = mol.GetNumHeavyAtoms()
            ring_bonds = sum([b.IsInRing() - b.GetIsAromatic()
                              for b in mol.GetBonds()])
            chiral_centers = len(Chem.FindMolChiralCenters(mol))

            scores.append(
                - 2.00 * np.power(total_atoms, 1.5)
                - 1.00 * np.power(ring_bonds, 1.5)
                - 2.00 * np.power(chiral_centers, 2.0)
            )

        sco = np.sum(scores) - 4.00 * np.power(necessary_reagent_atoms, 2.0)
        return sco / precursor['template_score']


    def reorder_precursors(self, precursors):
        scores = np.array([self.score_precursor(p) for p in precursors])
        indices = np.argsort(-scores)
        result = []
        for i, score in zip(indices, scores):
            result.append(precursors[i])
            result[-1]['score'] = score
        return result

    def get_priority(self, retroPrecursor, **kwargs):
        """Gets priority of given precursor based on heuristic and relevance.

        Args:
            retroPrecursor (RetroPrecursor): Precursor to calculate priority of.
            **kwargs: Unused.

        Returns:
            float: Priority score of precursor.
        """
        if not self._loaded:
            self.load_model()

        necessary_reagent_atoms = retroPrecursor.necessary_reagent.count('[') / 2.
        scores = []
        for smiles in retroPrecursor.smiles_list:
            # If buyable, basically free
            ppg = self.pricer.lookup_smiles(smiles, alreadyCanonical=True)
            if ppg:
                scores.append(- ppg / 1000.0)
                continue

            # Else, use heuristic
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

        sco = np.sum(scores) - 4.00 * np.power(necessary_reagent_atoms, 2.0)
        return sco / retroPrecursor.template_score

    def load_model(self):
        """Loads the Pricer used in the heuristic priority scoring."""
        self.pricer = Pricer()
        self.pricer.load()
        self._loaded = True
