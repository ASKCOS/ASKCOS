from makeit.prioritization.prioritizer import Prioritizer
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from makeit.utilities.io.logging import MyLogger
rank_context_prioritizer_loc = 'rank_context_prioritizer'


class RankContextPrioritizer(Prioritizer):

    def __init__(self):
        pass
    def get_priority(self, outcomes):
        '''
        Outcomes are the direct outcome of calling the single reaction evaluator
        '''
        return sorted(outcomes, key=lambda z: z['target']['rank'], reverse=True)

    def load_model(self):
        pass
