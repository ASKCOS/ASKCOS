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
        self.smiles_to_product = {}
        self.smiles_list_to_product = {}
        
    def add_product(self, product):
        '''
        Adds a product to the product set if it is a new product.
        Product type is ForwardProduct
        '''
        # Check if it is new or old
        try:
            index = self.smiles_to_product[product.smiles]
        except KeyError:
            try:
                index = self.smiles_list_to_product['.'.join(product.smiles_list)]
            except KeyError:
                #If neither has been encountered: add new product
                self.products.append(product)
                self.smiles_to_product[product.smiles] = len(self.products) - 1
                self.smiles_list_to_product['.'.join(product.smiles_list)] = len(self.products) - 1
                return
        
        self.products[index].template_ids += product.template_ids
        self.products[index].num_examples += product.num_examples
        
    def add_products(self, products):
        for product in products:
            self.add_product(product)
            
    def get_products(self):
        return self.products

class ForwardProduct:
    '''
    A class to store a single forward product for reaction enumeration
    '''
    def __init__(self, smiles_list = [], smiles = '', template_id = -1, num_examples = 0, 
                 edits = None, template_ids=None):
        self.smiles_list = smiles_list
        self.smiles = smiles
        self.template_ids = [template_id]
        if template_ids:
            self.template_ids = template_ids
        self.num_examples = num_examples
        self.edits = edits
    
    def get_edits(self):
        return self.edits

    def get_smiles(self):
        return self.smiles

    def as_dict(self):
        return {
            'smiles': self.smiles,
            'template_ids': [str(x) for x in self.template_ids],
            'num_examples': self.num_examples,
        }