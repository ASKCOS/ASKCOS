from makeit.prioritization.prioritizer import Prioritizer
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from makeit.utilities.io.logger import MyLogger
rank_context_prioritizer_loc = 'rank_context_prioritizer'


class RankContextPrioritizer(Prioritizer):
    """A context prioritizer that prioritizes on rank."""
    def __init__(self):
        """Initializes RankContextPrioritizer."""
        pass
    def get_priority(self, outcomes):
        """Gets priority of outcomes based on rank.

        Args:
            Outcomes (dict??): Direct outcome of calling the single reaction
                evaluator.

        Returns:
            float: Priority score (rank) of given outcomes.
        """
        return sorted(outcomes, key=lambda z: z['target']['rank'], reverse=True)

    def load_model(self):
        """Loads rank model.

        RankContextPrioritizer does not use a neural network, so this does
        nothing.
        """
        pass
