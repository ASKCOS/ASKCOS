import makeit.global_config as gc
import rdkit.Chem as Chem
from multiprocessing import Process, Manager, Queue
from makeit.utilities.io.logging import MyLogger
from makeit.utilities.io.model_loader import load_fastfilter, load_templatebased, load_templatefree
from makeit.utilities.parsing import parse_molecule_to_smiles, parse_list_to_smiles
from makeit.utilities.outcomes import summarize_reaction_outcome
from makeit.utilities.descriptors import edits_to_tensor
from celery.result import allow_join_result
from numpy import Inf

evaluator_loc = 'evaluator'


class Evaluator():
    '''
    The evaluator object uses any type of scorer with at least the 'get_scored_outcomes' method and any type of transformer
    with at least a 'get_outcomes' method.
    The outcomes are scored the results of the scoring can be returned in several ways.
    This is the major object that should be used in the synthetic evaluation direction



    reactants, target, context, transformer, scorer

            self.target_smiles = parse_molecule_to_smiles(target)
        self.reactants = Chem.MolFromSmiles(parse_list_to_smiles(reactants))
        self.reactants = clean_reactant_mapping(self.reactants)
        self.transformer = transformer
        self.scorer = scorer
        self.context = context
    '''

    def __init__(self, celery=False):

        self.celery = celery
        self.scorers = {}

    def evaluate(self, reactant_smiles, target, contexts, mincount=25, forward_scorer='', nproc=1, batch_size=250,
                 worker_no = 0, template_count=10000):
        with allow_join_result():
            target = Chem.MolToSmiles(Chem.MolFromSmiles(target))
            if not self.scorers:
                self.get_scorers(mincount, worker_no = worker_no)
            if not forward_scorer:
                MyLogger.print_and_log(
                    'Cannot evaluate a reaction without a forward scoring method. Exiting...', evaluator_loc, level=3)
            else:
                scorer = self.scorers[forward_scorer]

                all_outcomes = scorer.evaluate(reactant_smiles, contexts, batch_size=batch_size,
                                               template_prioritization=gc.popularity, nproc=nproc, soft_max=True,
                                               template_count=template_count)

                # output:
                # - top product for each context + score
                # - rank and score for target
                evaluation_results = []
                for i, outcomes in enumerate(all_outcomes):
                    evaluation_result = {'context': contexts[i]}
                    evaluation_result['top_product'] = {
                        'smiles': outcomes[0]['outcome'].smiles,
                        'template_ids': outcomes[0]['outcome'].template_ids,
                        'num_examples': outcomes[0]['outcome'].num_examples,
                        'rank': outcomes[0]['rank'],
                        'score': outcomes[0]['score'],
                        'prob': outcomes[0]['prob'],
                    }
                    evaluation_result['number_of_outcomes'] = len(outcomes)
                    found_target = False
                    for outcome in outcomes:
                        if target in outcome['outcome'].smiles:
                            found_target = True
                            evaluation_result['target'] = {
                                'smiles': outcome['outcome'].smiles,
                                'template_ids': outcome['outcome'].template_ids,
                                'num_examples': outcome['outcome'].num_examples,
                                'rank': outcome['rank'],
                                'score': outcome['score'],
                                'prob': outcome['prob'],
                            }
                    if not found_target:
                        MyLogger.print_and_log(
                            'Used template set did not find target in outcomes! Returning 0 as score.', evaluator_loc, level=2)
                        evaluation_result['target'] = {
                            'smiles': target,
                            'template_ids': [],
                            'num_examples': 0,
                            'rank': Inf,
                            'score': -Inf,
                            'prob': 0.0,
                        }
                    evaluation_results.append(evaluation_result)
            return evaluation_results

    def get_scorers(self, mincount, worker_no = 0):

        self.scorers[gc.fastfilter] = load_fastfilter(worker_no = worker_no)
        self.scorers[gc.templatebased] = load_templatebased(mincount=mincount, celery=self.celery, worker_no = worker_no)
        self.scorers[gc.templatefree] = load_templatefree(worker_no = worker_no)

if __name__ == '__main__':

    MyLogger.initialize_logFile()
    evaluator = Evaluator(celery=False)
    res = evaluator.evaluate('Cc1cccc(C)c1NC(=O)CCl', 'Cc1cccc(C)c1NC(=O)CN',
                             [[10.0, '', 'O=C(Cl)C(=O)Cl', '', 2.0, -1]], mincount=25, forward_scorer=gc.templatebased,
                             batch_size=1000, nproc=16)
    print res
    '''
    a = []
    for batch_size in range(250,1000,250):
        for nproc in range (4,8,2):
            a.append({'u':evaluator.evaluate('NC(=O)[C@H](CCC=O)N1C(=O)c2ccccc2C1=O', 'O=C1CC[C@H](N2C(=O)c3ccccc3C2=O)C(=O)N1', 
                             [(25, 'CN(C)C=O', '', '', 50, 60)], mincount = 25, forward_scorer=gc.templatebased,
                             batch_size = batch_size, nproc=nproc),
                      'b':batch_size,
                      'n':nproc})
    for aa in a:
        if aa['u'][0]['number_of_outcomes'] != 213:
            print('Incorrect outcomes for {} - {}'.format(aa['b'], aa['n']))
    '''
