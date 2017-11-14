import global_config as gc
import rdkit.Chem as Chem
from utilities.i_o.logging import MyLogger
from utilities.parsing import parse_molecule_to_smiles, parse_list_to_smiles
from utilities.outcomes import summarize_reaction_outcome
from utilities.descriptors import edits_to_tensor
from utilities.reactants import clean_reactant_mapping
template_evaluator_loc = 'evaluator'

class Evaluator():
    '''
    The evaluator object uses any type of scorer with at least the 'get_scored_outcomes' method and any type of transformer
    with at least a 'get_outcomes' method.
    The outcomes are scored the results of the scoring can be returned in several ways.
    This is the major object that should be used in the synthetic evaluation direction
    '''
    
    def __init__(self, reactants, target, context, transformer, scorer):
        self.target_smiles = parse_molecule_to_smiles(target)
        self.reactants = Chem.MolFromSmiles(parse_list_to_smiles(reactants))
        self.reactants = clean_reactant_mapping(self.reactants)
        self.transformer = transformer
        self.scorer = scorer
        self.scored_outcomes = []
        self.context = context
        
    def get_possible_products(self):
        return transformer.get_candidates(self.reactants, number_products)
    
    def run_scorer(self):
        forwardResults = self.transformer.get_outcomes(Chem.MolToSmiles(self.reactants))
        self.scored_outcomes = self.scorer.get_scored_outcomes(self.reactants, forwardResults, self.context, soft_max = True, sorted = True)
    
    def get_products_probability_cut(self, min_prob = 0.2):
        '''
        Get all products whose score exceeds a certain value.
        Outcomes must have been sorted
        Softmax has to have been applied 
        (both are default true in run_scorer)
        '''
        products = []
        if not self.scored_outcomes:
            self.run_scorer()
        index = 0
        (candidate_smiles, outcome, context, score) =  self.scored_outcomes[index]
          
        while score > min_prob:
            products.append((candidate_smiles, outcome, context, score))
            index += 1
            (candidate_smiles, outcome, context, score) =  self.scored_outcomes[index]
        
        return products
    
    def get_most_likely_product(self):
        return self.get_top_products(number = 1)
    
    def get_top_products(self, number = 3):
        '''
        Get the top-n products
        '''
        if self.scored_outcomes:
            return self.scored_outcomes[0:number]
        else:
            self.run_scorer()
            return self.scored_outcomes[0:number]
    
    def get_target_plausibility(self):
        target_score = 0.0
        target_found = False
        if self.scored_outcomes:
            print(self.target_smiles)
            for outcome in self.scored_outcomes:
                (candidate_smiles, outcome, context, score) = outcome
                if(self.target_smiles == candidate_smiles):
                    target_score += score
                    target_found = True
        else:
            self.run_scorer()
            for outcome in self.scored_outcomes:
                (candidate_smiles, outcome, context, score) = outcome
                if(self.target_smiles == candidate_smiles):
                    target_score += score
                    target_found = True
    
        if target_found:
            return target_score
        else:
            MyLogger.print_and_log('Target molecule not found among the options! Returning 0 as plausibility.', template_evaluator_loc, level = 1)
            return 0.0

            
if __name__ == '__main__':
    MyLogger.initialize_logFile()
    from parallelization.load_models import load_all
    from synthetic.forward_enumeration.forward_transformer import ForwardTransformer
    from scorer import Scorer
    ft = ForwardTransformer()
    ft.load()
    sc = Scorer()
    sc.load(folder = gc.PREDICTOR['trained_model_path'])
    eval = Evaluator("c1c(Cl)ccc(N)c1C(=O)c2ccccc2.NCC(=O)OCC","N1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13",('55','','','','',''), ft,sc)
    #eval = Evaluator("c1ccccc1C=O", "OC(c1ccccc1)C(=O)c2ccccc2", ('85',"O",'CCO','C#N','20','80'), ft, sc)
    #eval = Evaluator("CCCCC(=O)OCC.c1ccccc1O", "c1ccccc1OC(=O)CCCC", ('5',"O",'','','20','80'), ft, sc)
    print eval.get_top_products(number = 5)
    print eval.get_target_plausibility()
    print eval.get_products_probability_cut()