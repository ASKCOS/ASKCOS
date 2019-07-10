from makeit.prioritization.prioritizer import Prioritizer
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from makeit.utilities.io.logger import MyLogger
probability_context_prioritizer_loc = 'probability_context_prioritizer'


class ProbabilityContextPrioritizer(Prioritizer):
    """A context prioritizer that prioritizes on probability."""
    def __init__(self):
        """Initializes ProbabilityContextPrioritizer."""
        pass
    def get_priority(self, outcomes):
        """Gets priority of outcomes based on probability.

        Args:
            Outcomes (dict??): Direct outcome of calling the single reaction
                evaluator.

        Returns:
            float: Priority score (probability) of given outcomes.
        """
        return sorted(outcomes, key=lambda z: z['target']['prob'], reverse=True)

    def load_model(self):
        """Loads probability model.

        ProbabilityContextPrioritizer does not use a neural network, so this
        does nothing.
        """
        pass
