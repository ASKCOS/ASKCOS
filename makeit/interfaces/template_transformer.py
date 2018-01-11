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

        # Dictionary to keep track of ID -> index in self.templates
        self.id_to_index = {}
        
        return
    
    def get_prioritizers(self, prioritizers):
        '''
        Define which prioritization methods should be used. The prioritizers argument can contain 2 (retro - for precursor and templates) or 1 (synthetic - templates) element
        '''
        
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
    
    def get_prioritizers(self, args):
        '''
        Get the prioritization methods for the transformer (templates and/or precursors)
        '''
        raise NotImplementedError

    def top_templates(self):
        '''Generator to return only top templates. 
        Assumes templates are already sorted'''
        raise NotImplementedError
    
    def reorder(self):
        '''Reorder self.templates in descending popularity. Also builds id_tO_index table'''
        self.num_templates = len(self.templates)
        self.templates = sorted(self.templates, key=lambda z: z[
                                'count'], reverse=True)
        self.id_to_index = {template['_id']: i for i, template in enumerate(self.templates)}
        return

    def lookup_id(self, template_id):
        '''
        Find the reaction smarts for this template_id
        '''
        
        if not self.id_to_index: # need to build
            self.id_to_index = {template['_id']: i for (i, template) in enumerate(self.templates)}
        return self.templates[self.id_to_index[template_id]]
    
    #Define the methods that should be present in each transformer subclass.
    def load(self, chiral = False, lowe=False, refs=False, efgs=False, queue = None):
        raise NotImplementedError
    
    def get_outcomes(self, smiles, mincount, prioritizers, start_at = -1, end_at = -1, singleonly = False, stop_if = False, chiral=False):
        '''
        Performs a one-step retrosynthesis given a SMILES string of a
        target molecule by applying each transformation template
        sequentially.
        '''
        raise NotImplementedError
    
    def load_databases(self, chiral):
        
        raise NotImplementedError
    
    def apply_one_template(self, mol, smiles, template, singleonly = False, stop_if = False, chiral = False):
        '''
        Takes a mol object and applies a single template, returning
        a list of precursors or outcomes, depending on whether retro or 
        synthetic templates are used
        '''
        raise NotImplementedError
