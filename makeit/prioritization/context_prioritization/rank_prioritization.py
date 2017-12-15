from prioritization.prioritizer import Prioritizer
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from utilities.i_o.logging import MyLogger
rank_prioritizer_loc = 'probability_prioritizer'

class RankPrioritizer(Prioritizer):
    
    def __init__(self):
        MyLogger.print_and_log('Ranking contexts based on target rank.', rank_prioritizer_loc)
        
    def get_priority(self, outcomes):
        '''
        Outcomes are the direct outcome of calling the single reaction evaluator
        '''
        return sorted(outcomes, key = lambda z: z['target']['rank'], reverse = True)
    
    def load_model(self):
        pass