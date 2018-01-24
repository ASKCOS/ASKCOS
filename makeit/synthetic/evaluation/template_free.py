import makeit.global_config as gc
import rdkit.Chem as Chem
import json
import time
from pymongo import MongoClient
import numpy as np
import os
import sys
from celery.result import allow_join_result
from multiprocessing import Queue, Process, Manager
import Queue as VanillaQueue
from makeit.interfaces.scorer import Scorer

from makeit.synthetic.enumeration.results import ForwardResult, ForwardProduct
from makeit.synthetic.evaluation.template_based_aux import build
import makeit.utilities.contexts as context_cleaner


from makeit.utilities.parsing import parse_list_to_smiles
from makeit.utilities.io.logging import MyLogger
from makeit.utilities.reactants import clean_reactant_mapping
from askcos_site.askcos_celery.treeevaluator.forward_trans_worker import get_outcomes, template_count
from operator import itemgetter
template_nn_scorer_loc = 'template_based'


class TemplateFreeNeuralNetScorer(Scorer):

    def __init__(self, **kwargs):
        from makeit.synthetic.evaluation.rexgen_release.predict import TFFP
        self.model = TFFP() 

    def evaluate(self, reactants_smiles, contexts=None, **kwargs):
        '''Evaluation does not use context, but left as dummy pos var'''
        outcomes = self.model.predict(reactants_smiles, 
            top_n=kwargs.get('top_n', 1e10), num_core=kwargs.get('num_core', 8))
        if not outcomes:
            return [[{
                'rank': 1,
                'outcome': {'smiles': '', 'template_ids':[], 'num_examples':0},
                'score': 0,
                'prob': 0,
                }]]
         
        outcomes_to_ret = []
        for outcome in outcomes:
            outcomes_to_ret.append({ \
                'rank': outcome['rank'],
                'outcome': {
                    'smiles': max(outcome['smiles'].split('.'), key=len),
                    'smiles_full': outcome['smiles'],
                    'template_ids': [],
                    'num_examples': 0,
                },
                'score': outcome['score'],
                'prob': outcome['prob']})

        # Return in lists as if we received a list of contexts
        return [outcomes_to_ret]
                

if __name__ == '__main__':
    MyLogger.initialize_logFile()
    scorer = TemplateFreeNeuralNetScorer()
    res = scorer.evaluate('CN1C2CCC1CC(O)C2.O=C(O)C(CO)c1ccccc1')
    for re in res[0]:
        print re['outcome']['smiles'] + " {}".format(re['prob'])
    print 'done!'
