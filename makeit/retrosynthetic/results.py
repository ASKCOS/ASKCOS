import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import makeit.global_config as gc


class RetroResult:
    '''
    A class to store the results of a one-step retrosynthesis.
    '''

    def __init__(self, target_smiles):
        self.target_smiles = target_smiles
        self.precursors = []
        self.smiles_list_to_precursor = {}

    def add_precursor(self, precursor, prioritizer, **kwargs):
        '''
        Adds a precursor to the retrosynthesis result if it is a new and unique product
        '''
        try:
            index = self.smiles_list_to_precursor[
                '.'.join(precursor.smiles_list)]
        except KeyError:
            # If neither has been encountered: add new product
            precursor.prioritize(prioritizer, mode=kwargs.get('mode', gc.max))
            self.precursors.append(precursor)
            self.smiles_list_to_precursor[
                '.'.join(precursor.smiles_list)] = len(self.precursors) - 1
            return

        self.precursors[index].template_ids |= set(precursor.template_ids)
        self.precursors[index].num_examples += precursor.num_examples
        if self.precursors[index].template_score < precursor.template_score:
            self.precursors[index].template_score = precursor.template_score

    def return_top(self, n=50):
        '''
        Returns the top n precursors as a list of dictionaries, 
        sorted by descending score
        '''
        top = []
        for (i, precursor) in enumerate(sorted(self.precursors,
                                               key=lambda x: x.retroscore, reverse=True)):
            # Casts to float are necessary to maintain JSON serializability
            # when using celery
            top.append({
                'rank': i + 1,
                'smiles': '.'.join(precursor.smiles_list),
                'smiles_split': precursor.smiles_list,
                'score': float(precursor.retroscore),
                'num_examples': precursor.num_examples,
                'tforms': sorted(list(precursor.template_ids)),
                'template_score': float(precursor.template_score),
                'necessary_reagent': precursor.necessary_reagent,
            })
            if i + 1 == n:
                break
        return top


class RetroPrecursor:
    '''
    A class to store a single set of precursor(s) for a retrosynthesis
    does NOT contain the target molecule information
    '''

    def __init__(self, smiles_list=[], template_id=-1, template_score=1, num_examples=0, necessary_reagent=''):
        self.retroscore = 0
        self.num_examples = num_examples
        self.smiles_list = smiles_list
        self.template_ids = set([template_id])
        self.template_score = template_score
        self.necessary_reagent = necessary_reagent

    def prioritize(self, prioritizer, mode=gc.max):
        '''
        Calculate the score of this step as the worst of all precursors,
        plus some penalty for a large necessary_reagent
        '''
        self.retroscore = prioritizer.get_priority(self, mode=mode)
