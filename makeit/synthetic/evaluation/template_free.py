import makeit.global_config as gc
import rdkit.Chem as Chem
import json
import time
from pymongo import MongoClient
import numpy as np
import os
import sys
from makeit.interfaces.scorer import Scorer

from makeit.synthetic.enumeration.results import ForwardResult, ForwardProduct
from makeit.synthetic.evaluation.template_based_aux import build
import makeit.utilities.contexts as context_cleaner

from makeit.utilities.parsing import parse_list_to_smiles
from makeit.utilities.io.logging import MyLogger
from makeit.utilities.reactants import clean_reactant_mapping
from askcos_site.askcos_celery.treeevaluator.forward_trans_worker import get_outcomes, template_count
from operator import itemgetter
template_free_scorer_loc = 'template_based'


class TemplateFreeNeuralNetScorer(Scorer):

    def __init__(self, **kwargs):
        from makeit.synthetic.evaluation.rexgen_release.predict import TFFP
        self.model = TFFP() 

    def evaluate(self, reactants_smiles, contexts=[(20,'','','','','')], **kwargs):
        '''Evaluation does not use context, but left as dummy pos var'''

        all_outcomes = []
        for (T1, slvt1, rgt1, cat1, t1, y1) in contexts:
            this_reactants_smiles = reactants_smiles
            if slvt1:
                this_reactants_smiles += '.' + slvt1
            if rgt1:
                this_reactants_smiles += '.' + rgt1
            if cat1:
                this_reactants_smiles += '.' + cat1
            outcomes = self.model.predict(this_reactants_smiles, 
                top_n=kwargs.get('top_n', 1e10), num_core=kwargs.get('num_core', 8))
            if not outcomes:
                all_outcomes.append([{
                    'rank': 1,
                    'outcome': {'smiles': '', 'template_ids':[], 'num_examples':0},
                    'score': 0,
                    'prob': 0,
                    }])
                continue
             
            outcomes_to_ret = {}
            reactants_smiles_split = this_reactants_smiles.split('.')
            for outcome in outcomes:

                # Remove fragments that are unreacted
                smiles_list = set(outcome['smiles'].split('.'))
                for reactants_smiles_part in reactants_smiles_split:
                    if reactants_smiles_part in smiles_list:
                        smiles_list.remove(reactants_smiles_part)
                if not smiles_list:
                    continue # no reaction?

                # Canonicalize
                mol = Chem.MolFromSmiles('.'.join(smiles_list))
                if not mol:
                    MyLogger.print_and_log('Template free evaluator could not reparse product {}'.format('.'.join(smiles_list)), 
                        template_free_scorer_loc, 1)
                    continue
                smiles_list = Chem.MolToSmiles(mol).split('.')
                smiles = max(smiles_list, key=len) # NOTE: this is not great...byproducts may be longer

                if smiles in outcomes_to_ret:
                    outcomes_to_ret[smiles]['rank'] = min(outcomes_to_ret[smiles]['rank'], outcome['rank'])
                    outcomes_to_ret[smiles]['score'] = max(outcomes_to_ret[smiles]['score'], outcome['score'])
                    outcomes_to_ret[smiles]['prob'] += outcome['prob']
                else:
                    # Append outcome information
                    outcomes_to_ret[smiles] = {
                        'rank': outcome['rank'],
                        'outcome': {
                            'smiles': smiles,
                            # 'smiles_list': smiles_list,
                            'template_ids': [],
                            'num_examples': 0,
                        },
                        'score': float(outcome['score']),
                        'prob': float(outcome['prob'])
                    }

            # Renormalize and re-rank
            outcomes = sorted(outcomes_to_ret.values(), key=lambda x: x['prob'], reverse=True)
            total_prob = sum([outcome['prob'] for outcome in outcomes])
            for i, outcome in enumerate(outcomes):
                outcomes[i]['rank'] = i + 1
                outcomes[i]['prob'] = outcome['prob'] / total_prob

            all_outcomes.append(outcomes)
        # Return in lists as if we received a list of contexts
        return all_outcomes
                

if __name__ == '__main__':
    MyLogger.initialize_logFile()
    scorer = TemplateFreeNeuralNetScorer()
    res = scorer.evaluate('CN1C2CCC1CC(O)C2.O=C(O)C(CO)c1ccccc1')
    for re in res[0]:
        print(re['outcome']['smiles'] + " {}".format(re['prob']))
    print('done!')
