import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np

class ForwardResult:
    '''
    A class to store the results of a one-step forward synthesis.
    Should be used by any type of forward transformer/enumerator to store results and maintain uniformity!
    '''

    def __init__(self, smiles):
        self.smiles = smiles 
        self.products = []

    def add_product(self, product):
        '''
        Adds a product to the product set if it is a new product.
        Product type is ForwardProduct
        '''
        # Check if it is new or old
        for old_product in self.products:
            if product.smiles_list == old_product.smiles_list:
                # Just add this template_id and score
                old_product.template_ids |= set(product.template_ids)
                old_product.num_examples += product.num_examples
                return
        # New!
        self.products.append(product)

class ForwardProduct:
    '''
    A class to store a single forward product for reaction enumeration
    '''
    def __init__(self, smiles_list = [], template_id = -1, num_examples = 0):
        self.smiles_list = smiles_list
        self.template_ids = set([template_id])
        self.num_examples = num_examples