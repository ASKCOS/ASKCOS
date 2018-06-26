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

    def evaluate(self, reactant_smiles, target, contexts, forward_scorer='Template_Based',
                 worker_no=0, top_n=100, return_all_outcomes=False, chiral=False, **kwargs):
        with allow_join_result():
            target = Chem.MolToSmiles(Chem.MolFromSmiles(target), chiral)

            if not self.scorers:
                self.get_scorers(kwargs.get('mincount', 25), worker_no=worker_no)
            if not forward_scorer:
                MyLogger.print_and_log(
                    'Cannot evaluate a reaction without a forward scoring method. Exiting...', evaluator_loc, level=3)
            else:
                scorer = self.scorers[forward_scorer]
                if forward_scorer!= 'Fast_Filter':
                    all_outcomes = scorer.evaluate(reactant_smiles, contexts, **kwargs)
                else:
                    all_outcomes = scorer.evaluate(reactant_smiles, target, **kwargs)
                    # print(all_outcomes)
                    # wait = raw_input('press key to cont')
                # output:
                # - top product for each context + score
                # - rank and score for target
                evaluation_results = []
                for i, outcomes in enumerate(all_outcomes):
                    evaluation_result = {'context': contexts[i]}
                    if return_all_outcomes:
                        evaluation_result['outcomes'] = outcomes[:top_n]
                    
                    notadict = type(outcomes[0]['outcome']) != dict

                    evaluation_result['top_product'] = {
                        'smiles': outcomes[0]['outcome'].smiles if notadict else outcomes[0]['outcome']['smiles'],
                        'template_ids': outcomes[0]['outcome'].template_ids if notadict else outcomes[0]['outcome']['template_ids'],
                        'num_examples': outcomes[0]['outcome'].num_examples if notadict else outcomes[0]['outcome']['num_examples'],
                        'rank': outcomes[0]['rank'],
                        'score': outcomes[0]['score'],
                        'prob': outcomes[0]['prob'],
                    }
                    evaluation_result['number_of_outcomes'] = len(outcomes)
                    found_target = False
                    for outcome in outcomes:
                        outcome_smiles = outcome['outcome'].smiles if notadict else outcome['outcome']['smiles']
                        if target == outcome_smiles:
                            found_target = True
                            evaluation_result['target'] = {
                                'smiles': outcome['outcome'].smiles if notadict else outcome['outcome']['smiles'],
                                'template_ids': outcome['outcome'].template_ids if notadict else outcome['outcome']['template_ids'],
                                'num_examples': outcome['outcome'].num_examples if notadict else outcome['outcome']['num_examples'],
                                'rank': outcome['rank'],
                                'score': outcome['score'],
                                'prob': outcome['prob'],
                            }
                    if not found_target:
                        MyLogger.print_and_log(
                            'Target not found in outcomes! Returning 0 as score.', evaluator_loc, level=2)
                        MyLogger.print_and_log(
                            'The expected outcome of {} was {} ({})'.format(reactant_smiles, evaluation_result['top_product']['smiles'], 
                                evaluation_result['top_product']['prob']), evaluator_loc, level=2)
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

    def get_scorers(self, mincount, worker_no=0):

        self.scorers[gc.fastfilter] = load_fastfilter()
        self.scorers[gc.templatebased] = load_templatebased(
            mincount=mincount, celery=self.celery, worker_no=worker_no)
        self.scorers[gc.templatefree] = load_templatefree() # fast, one worker only

if __name__ == '__main__':

    MyLogger.initialize_logFile()
    evaluator = Evaluator(celery=False)
    # res = evaluator.evaluate('Cc1cccc(C)c1NC(=O)CCl.N', 'Cc1cccc(C)c1NC(=O)CN',
    #                          [[10.0, '', 'N', '', 2.0, -1]], mincount=25, forward_scorer=gc.templatebased,
    #                          batch_size=1000, nproc=16)
    # print res
    
    res = evaluator.evaluate('O=C1CCCCCCC1.OO', 'O=C1CCCCCCCO1', [(20,'','','','','')], forward_scorer=gc.templatefree)
    print(res)

    # ## the specificationo of other parameters should be automatically ignored if useing fastfilter
    # res = evaluator.evaluate('[ClH:8].[O:6]=[CH2:7].[s:1]1[cH:2][cH:3][cH:4][cH:5]1', 'O=Cc1cccs1',
    #                          [[10.0, '', 'O=C(Cl)C(=O)Cl', '', 2.0, -1]], mincount=25, forward_scorer=gc.fastfilter,
    #                          batch_size=1000, nproc=16)
    # print res
    
    # res = evaluator.evaluate('Cc1cccc(C)c1NC(=O)CCl', 'Cc1cccc(C)c1NC(=O)CN', ['unk.'], forward_scorer=gc.fastfilter)
    # print(res)
    
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
    c1cccc(C)c1NC(=O)CCl'''
