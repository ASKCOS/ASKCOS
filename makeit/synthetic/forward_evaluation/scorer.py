import global_config as gc
import rdkit.Chem as Chem   
import json
import numpy as np
import os
import sys
from forward_prediction_network import build
import utilities.contexts as context_cleaner
from utilities.outcomes import summarize_reaction_outcome
from utilities.descriptors import edits_to_vectors, edits_to_tensor
from utilities.parsing import parse_list_to_smiles
from utilities.i_o.logging import MyLogger
from operator import itemgetter
scorer_loc = 'scorer'

class Scorer:
    
    def __init__(self, done = None):
        self.model = None
        mol = Chem.MolFromSmiles('[C:1][C:2]')
        (a, _, b, _) = edits_to_vectors((['1'],[],[('1','2',1.0)],[]), mol)
        self.F_atom = len(a[0])
        self.F_bond = len(b[0])
        self.done = done
    
    def load(self, folder = ""):
        '''Load a neural network model'''
        if not folder:
            MyLogger.print_and_log('Cannot load neural network without the directory in which the parameters are saved. Exiting...', scorer_loc, level = 3)
        # Get model args
        ARGS_FPATH = os.path.join(folder, 'args.json')
        with open(ARGS_FPATH, 'r') as fid:
            args = json.load(fid)
    
        N_h2 = int(args['Nh2'])
        N_h1 = int(args['Nh1'])
        N_h3 = int(args['Nh3'])
        N_hf = int(args['Nhf'])
        l2v = float(args['l2'])
        lr = float(args['lr'])
        context_weight = float(args['context_weight'])
        enhancement_weight = float(args['enhancement_weight'])
        optimizer          = args['optimizer']
        inner_act          = args['inner_act']
        TARGET_YIELD       = False
        
        self.model = build(F_atom = self.F_atom, F_bond = self.F_bond, N_h1 = N_h1, 
                N_h2 = N_h2, N_h3 = N_h3, N_hf = N_hf, l2v = l2v, inner_act = inner_act,
                context_weight = context_weight, enhancement_weight = enhancement_weight, TARGET_YIELD = TARGET_YIELD,
                absolute_score = True)
    
        WEIGHTS_FPATH = os.path.join(folder, 'weights.h5')
        self.model.load_weights(WEIGHTS_FPATH, by_name = True)
        
        MyLogger.print_and_log('Scorer has been loaded.', scorer_loc)
        
        #multiprocessing notify done
        if self.done == None:
            pass
        else:
            self.done.value = 1
        
    def score_reaction(self, reactants, target_product, context):
        '''
        reactants: rdkit objects for the reactants (one object)
        target_product: rdkit object of the product (largest product: 1 molecule)
        context: full context
        '''
        prediction_context = context_cleaner.clean_context(context)
        prediction_context = context_cleaner.context_to_edit(prediction_context)
        prediction_edits = summarize_reaction_outcome(reactants, target_product)
        prediction_tensor = edits_to_tensor(prediction_edits, reactants, lenAtom = self.F_atom, lenBond = self.F_bond)
        full_tensor = prediction_tensor + prediction_context
        
        return self.model.predict(full_tensor)[0][0]
    
    def get_best_context_(self, reactants, target_product, contexts):
        '''
        From a list of contexts, get the context that results in the highest absolute score
        '''
        scored_contexts = []
        
        for context in contexts:
            score = self.score_reaction(reactants, target_product, context)

            scored_contexts.append((score, context))
            
        scored_contexts.sort(key=itemgetter(0), reverse=True)
        
        return scored_contexts[0]
    
    def get_scored_outcomes(self, reactants, forwardResults, context, soft_max = True, sorted = False):
        '''
        Of a list of possible forwardResults, score each possible outcome. Default uses softmax filter on scores.
        Each element in forwardResults should be of the type ForwardProduct. Reactants should be mapped.
        '''
        scored_outcomes = []
        for j, outcome in enumerate(forwardResults):
            outcome = Chem.MolFromSmiles(parse_list_to_smiles(outcome.smiles_list)) 
            try:
                           
                # Reduce to largest (longest) product only
                candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = gc.USE_STEREOCHEMISTRY)
                candidate_smiles = max(candidate_smiles.split('.'), key = len)
                outcome = Chem.MolFromSmiles(candidate_smiles)
                score = self.score_reaction(reactants, outcome, context)
                sys.stdout.flush()
                sys.stderr.flush()
                
                # Remove mapping before matching
                [a.ClearProp('molAtomMapNumber') for a in outcome.GetAtoms() if a.HasProp('molAtomMapNumber')] # remove atom mapping from outcome
                
                # Overwrite candidate_smiles without atom mapping numbers
                candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = gc.USE_STEREOCHEMISTRY)
                scored_outcomes.append((candidate_smiles, outcome, context, score))

            except Exception as e:
                [a.ClearProp('molAtomMapNumber') for a in outcome.GetAtoms() if a.HasProp('molAtomMapNumber')]
                candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = gc.USE_STEREOCHEMISTRY)
                MyLogger.print_and_log('Failed to score forwardResults for {}: {}.'.format(candidate_smiles, e), scorer_loc, level = 2)
                continue  
     
        
        if sorted:
            scored_outcomes.sort(key=itemgetter(3), reverse=True) 
            
        if soft_max:
            temp = scored_outcomes
            scored_outcomes = []
            scores = []
            info = []
            for (candidate_smiles, outcome, context, score) in temp:
                scores.append(score)
                info.append((candidate_smiles, outcome, context))
            scores = softmax(scores)
            if len(scores) == len(info):
                for i,score in enumerate(scores):
                    (candidate_smiles, outcome, context) = info[i]
                    scored_outcomes.append((candidate_smiles, outcome, context,score))
        
        return scored_outcomes

def softmax(x):
    """Compute softmax values for each sets of scores in x."""
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum()

            