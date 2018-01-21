from makeit.prioritization.prioritizer import Prioritizer
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from makeit.utilities.io.logging import MyLogger
probability_context_prioritizer_loc = 'probability_context_prioritizer'


class ProbabilityContextPrioritizer(Prioritizer):

    def __init__(self):
        pass
    def get_priority(self, outcomes):
        '''
        Outcomes are the direct outcome of calling the single reaction evaluator
        '''
        return sorted(outcomes, key=lambda z: z['target']['prob'], reverse=True)

    def load_model(self):
        pass
