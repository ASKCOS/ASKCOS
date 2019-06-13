import makeit.global_config as gc
import rdkit.Chem as Chem
import json
import time
from pymongo import MongoClient
import numpy as np
import os
import sys
from makeit.interfaces.scorer import Scorer
import makeit.utilities.contexts as context_cleaner
from makeit.utilities.io.logger import MyLogger
from operator import itemgetter
template_free_scorer_loc = 'template_free_scorer'


class TemplateFreeNeuralNetScorer(Scorer):
    """Template-free neural net evaluator.

    Attributes:
        model (TFFP): Template-free forward predictor.
    """
    def __init__(self, **kwargs):
        """Initializes TemplateFreeNeuralNetScorer.

        Args:
            **kwargs: Unused.
        """
        from makeit.synthetic.evaluation.rexgen_direct.predict import TFFP
        self.model = TFFP()

    def evaluate(self, reactants_smiles, contexts=[(20,'','','','','')], **kwargs):
        """Evaluates possible reaction outcomes for given reactants.

        Evaluation does not use context, but left as dummy pos var.

        Args:
            reactants_smiles (list of str??): SMILES string of reactants.
            contexts (list, optional): Unused.
                (default: {[(20,'','','','','')]})
            **kwargs: Additional optional parameters.

        Returns:
            list: Predicted reaction outcomes.
        """

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
                top_n=kwargs.get('top_n', 1e10))
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
                smiles_list = set(outcome['smiles'].split('.'))

                # Canonicalize
                smiles_canonical = set()
                for smi in smiles_list:
                    mol = Chem.MolFromSmiles(smi)
                    if not mol:
                        MyLogger.print_and_log('Template free evaluator could not reparse product {}'.format('.'.join(smiles_list)),
                        template_free_scorer_loc, 1)
                        continue
                    smiles_canonical.add(Chem.MolToSmiles(mol))


                # Remove unreacted frags
                smiles_canonical = smiles_canonical - set(reactants_smiles_split)
                if not smiles_canonical:
                    continue # no reaction?

                smiles = max(smiles_canonical, key=len) # NOTE: this is not great...byproducts may be longer
                if not smiles:
                    continue
                if smiles in outcomes_to_ret:
                    outcomes_to_ret[smiles]['rank'] = min(outcomes_to_ret[smiles]['rank'], outcome['rank'])
                    outcomes_to_ret[smiles]['score'] = np.log(np.exp(outcomes_to_ret[smiles]['score']) + np.exp(outcome['score']))
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
    import sys
    if len(sys.argv) > 1:
        react = str(sys.argv[1])
    else:
        react = 'CCCCO.CCCCBr'
    MyLogger.initialize_logFile()
    scorer = TemplateFreeNeuralNetScorer()
    res = scorer.evaluate(react)
    for re in res[0]:
        print(re['outcome']['smiles'] + " {}".format(re['prob']))
    print('done!')
