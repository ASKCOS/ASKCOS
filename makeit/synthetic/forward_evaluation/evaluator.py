import global_config as gc
import rdkit.Chem as Chem
from multiprocessing import Process, Manager, Queue
from utilities.i_o.logging import MyLogger
from utilities.i_o.model_loader import load_fastfilter, load_templatebased, load_templatefree
from utilities.parsing import parse_molecule_to_smiles, parse_list_to_smiles
from utilities.outcomes import summarize_reaction_outcome
from utilities.descriptors import edits_to_tensor
from celery.result import allow_join_result

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
    
    def __init__(self, celery = False):
        
        self.celery = celery
        self.scorers = {}
        
    def evaluate(self, reactant_smiles, target, contexts, mincount = 25, chiral = False,
                 forward_scorer='', nproc = 1, batch_size = 250):
        
        with allow_join_result():
            target = Chem.MolToSmiles(Chem.MolFromSmiles(target))     
            if not self.scorers:
                self.get_scorers(mincount, chiral)
            if not forward_scorer:
                MyLogger.print_and_log('Cannot evaluate a reaction without a forward scoring method. Exiting...', evaluator_loc, level = 3)
            else:
                scorer = self.scorers[forward_scorer]
                
                all_outcomes  =scorer.evaluate(reactant_smiles, contexts, batch_size = batch_size, 
                                template_prioritization=gc.popularity, nproc = nproc, soft_max = True)
                
                #output: 
                # - top product for each context + score
                # - rank and score for target
                evaluation_results = []
                for i,outcomes in enumerate(all_outcomes):
                    evaluation_result = {'context':contexts[i]}
                    evaluation_result['top_product']={
                                                      'smiles':outcomes[0]['outcome'].smiles,
                                                      'template_ids':outcomes[0]['outcome'].template_ids,
                                                      'num_examples':outcomes[0]['outcome'].num_examples,
                                                      'rank':outcomes[0]['rank'],
                                                      'score':str(outcomes[0]['score']),
                                                      'prob':str(outcomes[0]['prob']),
                                                      }
                    evaluation_result['number_of_outcomes']=len(outcomes)
                    for outcome in outcomes:
                        if target in outcome['outcome'].smiles:
                            evaluation_result['target']={
                                                        'smiles':outcome['outcome'].smiles,
                                                        'template_ids':outcome['outcome'].template_ids,
                                                        'num_examples':outcome['outcome'].num_examples,
                                                        'rank':outcome['rank'],
                                                        'score':str(outcome['score']),
                                                        'prob':str(outcome['prob']),
                                                        }
                    evaluation_results.append(evaluation_result)
            return evaluation_results           
                
            
    def get_scorers(self, mincount, chiral):
        
        self.scorers[gc.fastfilter] = load_fastfilter()
        self.scorers[gc.templatebased] = load_templatebased(mincount=mincount, chiral=chiral, celery = self.celery)
        self.scorers[gc.templatefree] = load_templatefree()
            
if __name__ == '__main__':
    MyLogger.initialize_logFile()
    evaluator = Evaluator(celery=False)
    print evaluator.evaluate('NC(=O)[C@H](CCC=O)N1C(=O)c2ccccc2C1=O', 'O=C1CC[C@H](N2C(=O)c3ccccc3C2=O)C(=O)N1', 
                             [(25, 'CN(C)C=O', '', '', 50, 60)], mincount = 25, forward_scorer=gc.templatebased,
                             batch_size = 50, nproc=1)