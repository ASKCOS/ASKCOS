import os
import global_config as gc
from i_o import arg_parser,name_parser, files
import rdkit.Chem as Chem
from i_o.logging import MyLogger
import i_o.model_loader as model_loader
from functions.treeBuilder import TreeBuilder
import sys
makeit_loc = 'makeit'

class MAKEIT:
    
    def __init__(self,TARGET, expansion_time, max_depth, max_branching, max_trees,retro_mincount,synth_mincount,
                 RANK_THRESHOLD_FOR_INCLUSION, PROB_THRESHOLD_FOR_INCLUSION,MAX_NUM_CONTEXTS,max_expansions,
                 max_ppg, min_trees_success, output_dir, nproc):
        
        self.TARGET = TARGET
        self.expansion_time = expansion_time
        self.max_depth = max_depth
        self.max_branching = max_branching
        self.max_trees = max_trees
        self.retro_mincount = retro_mincount
        self.synth_mincount = synth_mincount
        self.RANK_THRESHOLD_FOR_INCLUSION = RANK_THRESHOLD_FOR_INCLUSION
        self.PROB_THRESHOLD_FOR_INCLUSION = PROB_THRESHOLD_FOR_INCLUSION
        self.MAX_NUM_CONTEXTS = MAX_NUM_CONTEXTS
        self.max_expansions = max_expansions
        self.max_ppg = max_ppg
        self.min_trees_success = min_trees_success
        self.mol = name_parser.name_to_molecule(TARGET)
        self.smiles = Chem.MolToSmiles(self.mol)
        self.ROOT = files.make_directory(output_dir)
        self.case_dir = files.make_directory('{}/{}'.format(self.ROOT, self.TARGET))
        self.nproc = nproc
        
        
    def construct_synthesis_tree(self):
        treeBuilder = TreeBuilder(nb_workers = self.nproc, max_depth = self.max_depth, max_branching = self.max_branching, mincount = self.retro_mincount, 
                                  Pricer = self.models['pricer'], RetroTransformer = self.models['retro_transformer'],
                                  ContextModel = self.models['contextRecommender'], ForwardPredictor = self.models['forwardPredictor'], Scorer = self.models['scorer'])
        treeBuilder.build_tree(self.smiles)
        MyLogger.print_and_log(treeBuilder.tree_dict, makeit_loc)
        return treeBuilder
        
    
    def load_models(self):
        MyLogger.print_and_log('Loading models...', makeit_loc)
        databases = model_loader.load_Databases()
        retroTransformer = model_loader.load_Retro_Transformer(databases['Retro_Database'], mincount_retro = self.retro_mincount)
        synthTransformer = model_loader.load_Forward_Transformer(databases['Synth_Database'], mincount_synth = self.synth_mincount)
        pricer = model_loader.load_Pricer(databases['Chemical_Database'], databases['Buyable_Database'])
        forwardPredictor = model_loader.load_Forward_Predictor(databases['Synth_Database'], mincount_synth = self.synth_mincount)
        scorer = model_loader.load_Scorer()
        contextRecommender = model_loader.load_Context_Recommender(databases['Reaction_Database'], databases['Instance_Database'], 
                                                                   databases['Chemical_Database'], databases['Solvent_Database'], max_total_contexts = self.MAX_NUM_CONTEXTS)
        self.models ={
            'retro_transformer':retroTransformer,
            'synthetic_transformer':synthTransformer,
            'pricer':pricer,
            'forwardPredictor':forwardPredictor,
            'contextRecommender':contextRecommender,
            'scorer':scorer
            }
        MyLogger.print_and_log('All models loaded.', makeit_loc)
        
def find_synthesis():
    args = arg_parser.get_args()
    makeit = MAKEIT(args.TARGET, args.expansion_time, args.max_depth, args.max_branching,
                    args.max_trees,args.retro_mincount,args.synth_mincount, args.rank_threshold, 
                    args.p_threshold, args.max_contexts, args.max_expansions, args.max_ppg,
                    args.min_trees_success, args.output, args.nproc)
    
    MyLogger.initialize_logFile(makeit.ROOT, makeit.case_dir)
    
    makeit.load_models()
    
    tree = makeit.construct_synthesis_tree()
    

if __name__ == '__main__':
    find_synthesis()