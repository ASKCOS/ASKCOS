import os
import time
import global_config as gc
from utilities.i_o import arg_parser,name_parser, files
import rdkit.Chem as Chem
from utilities.i_o.logging import MyLogger
from retro_synthetic.tree_builder import TreeBuilder
from askcos_site.askcos_celery.treebuilder.tb_coordinator import get_buyable_paths
from askcos_site.askcos_celery.treeevaluator.tree_evaluation_coordinator import evaluate_trees
from synthetic.forward_evaluation.tree_evaluator import TreeEvaluator
import sys
makeit_loc = 'makeit'

class MAKEIT:
    '''
    Main application for running the make-it program. 
    Proposes potential synthetic routes to a desired target compound in two steps:
     - Building a retro synthetic tree and extracting buyable routes
     - Evaluation the likelihood of succes of each of the reactions in the found buyable routes
     - Returns all (or some) of the likely synthetic routes
    '''
    def __init__(self,TARGET, expansion_time, max_depth, max_branching, max_trees,retro_mincount, retro_mincount_c, 
                 synth_mincount, rank_threshold_inclusion, prob_threshold_inclusion,max_total_contexts,
                 max_ppg, output_dir, chiral, nproc, celery, context_recommender, forward_scoring_method,
                 tree_scoring_method, context_prioritization, template_prioritization, precursor_prioritization):
        
        self.TARGET = TARGET
        self.expansion_time = expansion_time
        self.max_depth = max_depth
        self.max_branching = max_branching
        self.max_trees = max_trees
        self.context_recommender = context_recommender
        self.forward_scoring_method = forward_scoring_method
        self.tree_scoring_method = tree_scoring_method
        self.context_prioritization = context_prioritization
        self.template_prioritization = template_prioritization
        self.precursor_prioritization = precursor_prioritization
        self.retro_mincount = retro_mincount
        self.retro_mincount_c = retro_mincount_c
        self.synth_mincount = synth_mincount
        self.rank_threshold_inclusion = rank_threshold_inclusion
        self.prob_threshold_inclusion = prob_threshold_inclusion
        self.max_total_contexts = max_total_contexts
        self.max_ppg = max_ppg
        self.mol = name_parser.name_to_molecule(TARGET)
        self.smiles = Chem.MolToSmiles(self.mol)
        self.ROOT = files.make_directory(output_dir)
        self.case_dir = files.make_directory('{}/{}'.format(self.ROOT, self.TARGET))
        self.nproc = nproc
        self.celery = celery
        self.chiral = chiral
        self.known_bad_reactions = []
        
    def construct_buyable_trees(self):
        
        if self.celery:#Call celery worker
            working = time.time()
            res = get_buyable_paths.apply_async(args = (self.smiles, self.template_prioritization, self.precursor_prioritization),
                                                          kwargs = {'mincount':self.retro_mincount, 'max_branching':self.max_branching,
                                                                    'max_depth':self.max_depth, 'max_ppg': self.max_ppg, 'max_time':self.expansion_time,
                                                                    'max_trees':self.max_trees, 'known_bad_reactions':self.known_bad_reactions,
                                                                    'chiral':self.chiral})
            
            while not res.ready():
                if int(time.time() - working)%10 == 0:
                    MyLogger.print_and_log('Building trees...', makeit_loc)
                time.sleep(1)
            buyable_trees = res.get()
        else:#Create tree builder object and run it
            treeBuilder = TreeBuilder(celery = self.celery, mincount = self.retro_mincount, mincount_c = self.retro_mincount_c, chiral = self.chiral)
        
            buyable_trees = treeBuilder.get_buyable_paths(self.smiles, template_prioritization = self.template_prioritization,
                                                          precursor_prioritization = self.precursor_prioritization, nproc = self.nproc,
                                                          max_depth=self.max_depth, max_branching=self.max_branching, max_ppg = self.max_ppg,
                                                          mincount = self.retro_mincount, chiral = self.chiral, max_trees=self.max_trees,
                                                          known_bad_reactions=self.known_bad_reactions, expansion_time = self.expansion_time)
        
        return buyable_trees
    
    def evaluate_synthesis_trees(self, trees):
        if self.celery:#Call celery worker
            working = time.time()
            res = evaluate_trees.apply_async(args= (trees,), kwargs = {'context_scoring_method':self.context_prioritization,
                                                                            'context_recommender':self.context_recommender,
                                                                            'forward_scoring_method':self.forward_scoring_method,
                                                                            'tree_scoring_method':self.tree_scoring_method,
                                                                            'rank_threshold':self.rank_threshold_inclusion,
                                                                            'prob_threshold':self.prob_threshold_inclusion,
                                                                            'mincount':self.synth_mincount, 
                                                                            'batch_size':500, 'n':self.max_total_contexts})
            while not res.ready():
                if int(time.time() - working)%10 == 0:
                    MyLogger.print_and_log('Evaluating trees...', makeit_loc)
                time.sleep(1)
            evaluated_trees = res.get()
        else:#Create a tree evaluation object and run it
            if self.forward_scoring_method == gc.templatebased:
                #nproc = number of parallel forward enumeration workers
                #nproc_t = number of trees to be evaluated in parallel.
                #Only use an nproc different from 1 if using the template base forward evaluation method. Otherwise
                #evaluation is fast enough to do without additional parallelization
                if len(trees)>2:
                    parallel = True
                    nproc_t = 2
                    nproc = self.nproc/2
                else:
                    nproc_t = 1
                    nproc = self.nproc
                    parallel = False
            else:
                nproc_t = self.nproc
                nproc = 1
            treeEvaluator = TreeEvaluator(celery = False, context_recommender = self.context_recommender)
            evaluated_trees = treeEvaluator.evaluate_trees(trees, context_recommender = self.context_recommender, context_scoring_method=self.context_prioritization,
                                                    forward_scoring_method =self.forward_scoring_method, tree_scoring_method=self.tree_scoring_method,
                                                    rank_threshold =self.rank_threshold_inclusion, prob_threshold= self.prob_threshold_inclusion,
                                                    mincount = self.synth_mincount, batch_size = 500, n = self.max_total_contexts, nproc_t = nproc_t,
                                                    nproc = nproc, parallel = parallel) 
        plausible_trees = []
        print evaluated_trees
        for tree in evaluated_trees:
            if tree['plausible']:
                plausible_trees.append(tree)
                
        if plausible_trees:
            MyLogger.print_and_log('Feasible synthesis route discovered!', makeit_loc)
        else:
            MyLogger.print_and_log('No feasible routes from buyables have been discovered. Consider changing inclusion thesholds.', makeit_loc)
        
        return plausible_trees
    
def find_synthesis():
   
    args = arg_parser.get_args()
    makeit = MAKEIT(args.TARGET, args.expansion_time, args.max_depth, args.max_branching,
                    args.max_trees,args.retro_mincount,args.retro_mincount_c,args.synth_mincount,
                    args.rank_threshold,args.prob_threshold, args.max_contexts, args.max_ppg,
                    args.output, args.chiral, args.nproc, args.celery, args.context_recommender,
                    args.forward_scoring, args.tree_scoring, args.context_prioritization,
                    args.template_prioritization, args.precursor_prioritization)
    MyLogger.initialize_logFile(makeit.ROOT, makeit.case_dir)
    #only load modules if not using celery! Modules are preloaded in celery.

    trees = makeit.construct_buyable_trees()
    MyLogger.print_and_log('MAKEIT generated {} buyable tree(s) that meet(s) all constraints.'.format(len(trees)), makeit_loc)
    feasible_trees = makeit.evaluate_synthesis_trees(trees)
    MyLogger.print_and_log('MAKEIT found {} tree(s) that are(is) likely to result in a successful synthesis.'.format(len(feasible_trees)),makeit_loc)
    print feasible_trees
    
if __name__ == '__main__':
    find_synthesis()