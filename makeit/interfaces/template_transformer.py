from __future__ import print_function

class TemplateTransformer(object):
    '''
    The Transformer class defines an object which can be used to perform
    one-step retrosyntheses for a given molecule.
    '''

    def __init__(self, mincount=4, TEMPLATE_DB = None, loc = False, done = None):
        '''
        Initialize a transformer.
        TEMPLATE_DB: indicate the database you want to use (def. none)
        loc: indicate that local file data should be read instead of online data (def. false)
        '''
        
        raise NotImplementedError
        
    def dump_to_file(self, file_name):
        '''
        Write the template database to a file, of which the path in specified in the general configuration
        '''
        raise NotImplementedError
    
    def load_from_file(self, file_name):
        '''
        Read the template database from a previously saved file, of which the path is specified in the general
        configuration
        '''
        raise NotImplementedError
    
    def reorder(self):
        '''
        Re-orders the list of templates (self.templates) according to 
        field 'count' in descending order. This means we will apply the
        most popular templates first
        '''
        raise NotImplementedError

    def top_templates(self):
        '''Generator to return only top templates. 
        Assumes templates are already sorted'''
        raise NotImplementedError
    
    def lookup_id(self, template_id):
        '''
        Find the reaction smarts for this template_id
        '''
        raise NotImplementedError
    
    #Define the methods that should be present in each transformer subclass.
    def load(self, lowe=False, refs=False, efgs=False, queue = None):
        raise NotImplementedError
    
    def perform_transformations(self, smiles, singleonly = False, stop_if = False):
        '''
        Performs a one-step retrosynthesis given a SMILES string of a
        target molecule by applying each transformation template
        sequentially.
        '''
        raise NotImplementedError
    
    def load_databases(self):
        
        raise NotImplementedError
    
    def apply_one_template(self, mol, smiles, template, singleonly = False, stop_if = False, chiral = False):
        '''
        Takes a mol object and applies a single template, returning
        a list of precursors or outcomes, depending on whether retro or 
        synthetic templates are used
        '''
        raise NotImplementedError