from prioritization.prioritizer import Prioritizer
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from utilities.i_o.logging import MyLogger
probability_prioritizer_loc = 'probability_prioritizer'

class ProbabilityPrioritizer(Prioritizer):
    
    def __init__(self):
        MyLogger.print_and_log('Ranking contexts based on target probability.', probability_prioritizer_loc)
        
    def get_priority(self, outcomes):
        '''
        Outcomes are the direct outcome of calling the single reaction evaluator
        '''
        return sorted(outcomes, key = lambda z: z['target']['prob'], reverse = True)
        
    def load_model(self):
        pass