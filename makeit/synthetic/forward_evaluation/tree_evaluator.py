import global_config as gc
from multiprocessing import Process, Manager, Queue
from evaluator import Evaluator
import Queue as VanillaQueue
import time
from utilities.i_o.logging import MyLogger
from askcos_site.askcos_celery.contextrecommender.cr_coordinator import get_context_recommendations
from askcos_site.askcos_celery.treeevaluator.scoring_coordinator import evaluate
from prioritization.context_prioritization.probability_prioritization import ProbabilityPrioritizer
from prioritization.context_prioritization.rank_prioritization import RankPrioritizer
from prioritization.default_prioritizer import DefaultPrioritizer
treeEvaluator_loc = 'tree_evaluator'

class TreeEvaluator():
    '''
    Class for the evaluation of the found retrosynthetic tree
    '''
    def __init__(self, celery = False, forward_transformer = None, scorer = None, context_recommender = None, 
                 tree_dict = None, rank_inclusion = 10, prob_inclusion = 0.2, max_contexts = 10, single_solv = True,
                 with_smiles = True, recommender = None, batch_size = 500):
        self.celery = celery
        self.single_solv = single_solv
        self.with_smiles = with_smiles
        self.batch_size = batch_size
        self.evaluator = Evaluator(celery = self.celery)
        self.reset()
        self.rank_threshold = rank_inclusion
        self.prob_threshold = prob_inclusion
        
        if not self.celery:
            if context_recommender == None: #and not fast filter
                MyLogger.print_and_log('Context recommendation method required to evaluate the tree.', treeEvaluator_loc, level = 3)
        
        if self.celery:
            def get_contexts(rxn, n):
                res = get_context_recommendations.apply_async(args = (rxn), 
                                                              kwargs= {'n':n, 'singleSlvt':self.single_solv, 
                                                                       'with_smiles':self.with_smiles, 
                                                                       'recommender':self.recommender})
                while not res.ready():
                    time.sleep(0.2)
                return res.get()
        else:
            def get_contexts(rxn, n):
                return context_recommender.get_n_conditions(rxn, n=n, singleSlvt=self.single_solv, with_smiles=self.with_smiles)
        
        self.get_contexts = get_contexts
        
        if self.celery:
            def evaluate_reaction(reactant_smiles, contexts):
                res = evaluate.apply_async(args = (reactant_smiles, target, contexts), 
                                           kwargs = {'mincount': self.mincount, 'chiral' : self.chiral, 
                                                     'forward_scorer': self.forward_scorer})
                while not res.ready():
                    time.sleep(0.2)
                return res.get()
        else:
            def evaluate_reaction(reactant_smiles, target, contexts):
                return self.evaluator.evaluate(reactant_smiles, target, contexts, mincount = self.mincount, chiral = self.chiral,
                                               forward_scorer=self.forward_scorer, nproc = self.nproc, batch_size = self.batch_size)
        
        self.evaluate_reaction = evaluate_reaction
        
    def get_context_prioritizer(self, context_method):
        if context_method == gc.probability:
            self.context_prioritizer = ProbabilityPrioritizer()
        elif context_method == gc.rank:
            self.context_prioritizer = RankPrioritizer()
        else:
            MyLogger.print_and_log('Specified prioritization method does not exist. Using default method.', treeEvaluator_loc, level =1)
            self.context_prioritizer = DefaultPrioritizer()
        self.context_prioritizer.load_model()
        
        
    def get_top_context(self, evaluation):
        return self.context_prioritizer.get_priority(evaluation)[0]
    
    def reset(self):
       self.evaluation_dict = {}
       
    def is_plausible(self, result):
        prob = result['target']['prob']
        rank = result['target']['rank']
        return prob > self.prob_threshold and rank < self.rank_threshold
    
    def evaluate_tree(self, tree, context_method = '', scoring_method = '', rank_threshold = 5, prob_threshold = 0.2, is_target = False, 
                      mincount = 25, nproc = 1,batch_size = 500, chiral = False, n = 10):
        if is_target:
            self.reset()
            self.get_context_prioritizer(context_method)
            self.rank_threshold = rank_threshold
            self.prob_threshold = prob_threshold
            self.mincount = mincount
            self.nproc = nproc
            self.batch_size = batch_size
            self.chiral = chiral
            self.forward_scorer = scoring_method
        
        if not tree['children']:
            #Reached the end of the synthesis tree -> Stop
            return 1.0
        else:
            target = tree['smiles']
            reaction = tree['children'][0]
            reactants = [child['smiles'] for child in reaction['children']]
            reaction_smiles = reaction['smiles']
            necessary_reagent = reaction['necessary_reagent']
            ###############################################################
            #If reaction encountered before: get data from dict.
            ###############################################################
            if reaction_smiles in self.evaluation_dict:
                evaluation = self.evaluation_dict[reaction_smiles]
            ###############################################################
            #Otherwise create data
            ###############################################################
            else:
                if necessary_reagent:
                    contexts = self.get_contexts(reaction_smiles, 1)
                else:
                    contexts = self.get_contexts(reaction_smiles, n)
                evaluation = self.evaluate_reaction('.'.join(reactants), target, contexts)
                self.evaluation_dict[reaction_smiles] = evaluation
            ###############################################################
            #Process data
            ###############################################################
            
            if len(evaluation) == 1:
                top_result == evaluation[0]
            else:
                top_result = self.get_top_context(evaluation)
            tree['evaluation'] = top_result
            MyLogger.print_and_log('Evaluated reaction: {} - ranked {} with a {}% probability.'\
                                    .format(reaction_smiles, top_result['target']['rank'], float(top_result['target']['prob'])*100.0),\
                                    treeEvaluator_loc)
            plausible = self.is_plausible(top_result)
            all_children_plausible = True
            for child in reaction['children']:
                child_plausible =self.evaluate_tree(child)
        
                if not child_plausible:
                    all_children_plausible = False
                    
            if all_children_plausible and is_target:
                MyLogger.print_and_log('Found a fully plausible tree!', treeEvaluator_loc)
            elif is_target:
                MyLogger.print_and_log('Evaluated tree has unfeasible children.', treeEvaluator_loc)
            
            
            
            
            return plausible
    
if __name__ == '__main__':
    from synthetic.context.nn_context_recommender import NNContextRecommender
    cr = NNContextRecommender()
    cr.load_nn_model(model_path = gc.CONTEXT_REC['model_path'], info_path = gc.CONTEXT_REC['info_path'])
    ev = TreeEvaluator(context_recommender = cr)
    tree = {'is_chemical': True, 'smiles': 'CCC(=O)OCCCNCC', 'ppg': 0.0, 'id': 1, 'children': [{'info': '', 'smiles': 'CCC(=O)O.CCNCCCO>>CCC(=O)OCCCNCC', 'is_reaction': True, 'id': 2, 'num_examples': 26551, 'children': [{'is_chemical': True, 'smiles': 'CCC(=O)O', 'ppg': 1.0, 'id': 3, 'children': []}, {'is_chemical': True, 'smiles': 'CCNCCCO', 'ppg': 0.0, 'id': 4, 'children': [{'info': '', 'smiles': 'CC=O.NCCCO>>CCNCCCO', 'is_reaction': True, 'id': 381, 'num_examples': 3977, 'children': [{'is_chemical': True, 'smiles': 'CC=O', 'ppg': 1.0, 'id': 38, 'children': []}, {'is_chemical': True, 'smiles': 'NCCCO', 'ppg': 1.0, 'id': 320, 'children': []}], 'necessary_reagent': u''}]}], 'necessary_reagent': u''}]}
    res = ev.evaluate_tree(tree, gc.probability, gc.templatebased, is_target=True)