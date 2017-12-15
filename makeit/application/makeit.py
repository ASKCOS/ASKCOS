import os
import time
import global_config as gc
from utilities.i_o import arg_parser,name_parser, files
import rdkit.Chem as Chem
from utilities.i_o.logging import MyLogger
import utilities.i_o.model_loader as model_loader
from retro_synthetic.tree_builder import TreeBuilder
from synthetic.forward_evaluation.tree_evaluator import TreeEvaluator
import sys
makeit_loc = 'makeit'

class MAKEIT:
    
    def __init__(self,TARGET, expansion_time, max_depth, max_branching, max_trees,retro_mincount,synth_mincount,
                 rank_threshold_inclusion, prob_threshold_inclusion,max_total_contexts,max_expansions,
                 max_ppg, min_trees_success, output_dir, chiral, nproc, celery):
        
        self.TARGET = TARGET
        self.expansion_time = expansion_time
        self.max_depth = max_depth
        self.max_branching = max_branching
        self.max_trees = max_trees
        self.retro_mincount = retro_mincount
        self.synth_mincount = synth_mincount
        self.rank_threshold_inclusion = rank_threshold_inclusion
        self.prob_threshold_inclusion = prob_threshold_inclusion
        self.max_total_contexts = max_total_contexts
        self.max_expansions = max_expansions
        self.max_ppg = max_ppg
        self.min_trees_success = min_trees_success
        self.mol = name_parser.name_to_molecule(TARGET)
        self.smiles = Chem.MolToSmiles(self.mol)
        self.ROOT = files.make_directory(output_dir)
        self.case_dir = files.make_directory('{}/{}'.format(self.ROOT, self.TARGET))
        self.nproc = nproc
        self.celery = celery
        self.chiral = chiral
        
        
    def construct_synthesis_tree(self):
        
        treeBuilder = TreeBuilder(retroTransformer = self.models['retro_transformer'], pricer = self.models['pricer'], 
                                  max_branching = self.max_branching, max_depth = self.max_depth, expansion_time = self.expansion_time,
                                  celery = self.celery)
        
        buyable_trees = treeBuilder.get_buyable_paths(self.smiles)
        
        return buyable_trees
    
    def evaluate_synthesis_tree(self, tree):
        
        treeEvaluator = TreeEvaluator(forward_transformer = self.models['synthetic_transformer'], scorer = self.models['scorer'], context_recommender = self.models['context_recommender'],
                                        tree_dict = tree, rank_inclusion = self.rank_threshold_inclusion, prob_inclusion = self.prob_threshold_inclusion, max_contexts = self.max_total_contexts)
        
        feasible = treeEvaluator.evaluate_tree() 
        
        if feasible:
            MyLogger.print_and_log('Feasible synthesis route discovered!', makeit_loc)
        else:
            MyLogger.print_and_log('No feasible routes from buyables have been discovered. Consider changing inclusion thesholds.', makeit_loc)
        
        MyLogger.print_and_log(treeEvaluator.plausible_tree, makeit_loc)
        return treeEvaluator.plausible_tree
    
    def load_modules(self):
        self.models = model_loader.load_all(retro_mincount = self.retro_mincount, synth_mincount = self.synth_mincount, max_contexts = self.max_total_contexts)
        
def find_synthesis():
    args = arg_parser.get_args()
    makeit = MAKEIT(args.TARGET, args.expansion_time, args.max_depth, args.max_branching,
                    args.max_trees,args.retro_mincount,args.synth_mincount, args.rank_threshold, 
                    args.p_threshold, args.max_contexts, args.max_expansions, args.max_ppg,
                    args.min_trees_success, args.output, args.chiral, args.nproc, args.celery)
    
    MyLogger.initialize_logFile(makeit.ROOT, makeit.case_dir)
    
    #only load modules if not using celery! Modules are preloaded in celery.
    if not args.celery:
        makeit.load_modules()
    
    tree = makeit.construct_synthesis_tree()
    plausbible_tree = makeit.evaluate_synthesis_tree(tree)

if __name__ == '__main__':
    find_synthesis()