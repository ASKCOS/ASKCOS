from __future__ import print_function
import global_config as gc
from global_config import USE_STEREOCHEMISTRY
import rdkit.Chem as Chem          
from rdkit.Chem import AllChem
import numpy as np
from functools import partial # used for passing args to multiprocessing
from utilities.i_o.logging import MyLogger
from utilities.reactants import clean_reactant_mapping
from retro_enumeration import *
from pymongo import MongoClient
from interfaces.template_transformer import TemplateTransformer
from rdchiral.main import rdchiralRun
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
retro_transformer_loc = 'retro_transformer'

class RetroTransformer(TemplateTransformer):
    '''
    The Transformer class defines an object which can be used to perform
    one-step retrosyntheses for a given molecule.
    '''

    def __init__(self, prioritizer = None, mincount=4, TEMPLATE_DB = None, loc = False, done = None):
        '''
        Initialize a transformer.
        TEMPLATE_DB: indicate the database you want to use (def. none)
        loc: indicate that local file data should be read instead of online data (def. false)
        '''
        self.TEMPLATE_DB = TEMPLATE_DB
        self.done = done
        self.mincount = mincount
        self.templates = []
        self.id_to_index = {}
        if not prioritizer:
            MyLogger.print_and_log('Cannot run the Retro-Transformer without a prioritization method. Exiting...', retro_transformer_loc, level = 3)
        self.prioritizer = prioritizer
        
    def load(self, lowe=False, refs=False, efgs=False):
        '''
        Loads and parses the template database to a useable one
        '''
        # Save collection TEMPLATE_DB
        if not self.TEMPLATE_DB:
            self.load_databases()
        
        if self.mincount and 'count' in self.TEMPLATE_DB.find_one(): 
            filter_dict = {'count': { '$gte': self.mincount}}
        else: 
            filter_dict = {}

        # Look for all templates in collection
        to_retrieve = ['_id', 'reaction_smarts', 'necessary_reagent', 'count', 'intra_only']
        if refs:
            to_retrieve.append('references')
        if efgs:
            to_retrieve.append('efgs')
        for document in self.TEMPLATE_DB.find(filter_dict, to_retrieve):
            # Skip if no reaction SMARTS
            if 'reaction_smarts' not in document: continue
            reaction_smarts = str(document['reaction_smarts'])
            if not reaction_smarts: continue

            # Define dictionary
            template = {
                'name':                 document['name'] if 'name' in document else '',
                'reaction_smarts':      reaction_smarts,
                'incompatible_groups':  document['incompatible_groups'] if 'incompatible_groups' in document else [],
                'reference':            document['reference'] if 'reference' in document else '',
                'references':           document['references'] if 'references' in document else [],
                'rxn_example':          document['rxn_example'] if 'rxn_example' in document else '',
                'explicit_H':           document['explicit_H'] if 'explicit_H' in document else False,
                '_id':                  document['_id'] if '_id' in document else -1,
                'product_smiles':       document['product_smiles'] if 'product_smiles' in document else [], 
                'necessary_reagent':    document['necessary_reagent'] if 'necessary_reagent' in document else '',       
                'efgs':                 document['efgs'] if 'efgs' in document else None,
                'intra_only':           document['intra_only'] if 'intra_only' in document else False,
            }

            # Frequency/popularity score
            if 'count' in document: 
                template['count'] = document['count']
            elif 'popularity' in document:
                template['count'] = document['popularity']
            else:
                template['count'] = 1

            # Define reaction in RDKit and validate
            try:
                # Force reactants and products to be one molecule (not really, but for bookkeeping)
                reaction_smarts_retro = '(' + reaction_smarts.replace('>>', ')>>(') + ')'
                rxn = AllChem.ReactionFromSmarts(str(reaction_smarts_retro))
                #if rxn.Validate() == (0, 0):
                if rxn.Validate()[1] == 0: 
                    template['rxn'] = rxn
                else:
                    template['rxn'] = None
            except Exception as e:
                MyLogger.print_and_log('Couldnt load retro: {}: {}'.format(reaction_smarts_retro, e),retro_transformer_loc, level=1)
                template['rxn'] = None
            
            if not template['rxn']: continue
            # Add to list
            self.templates.append(template)

        self.num_templates = len(self.templates)
        self.reorder()
        
        
        MyLogger.print_and_log('Retro-synthetic transformer has been loaded.', retro_transformer_loc)
        
        #multiprocessing notify done
        if self.done == None:
            pass
        else:
            self.done.value = 1

    def perform_transformations(self, smiles, singleonly = False, stop_if = False):
        '''
        Performs a one-step retrosynthesis given a SMILES string of a
        target molecule by applying each transformation template
        sequentially.
        '''

        # Define mol to operate on
        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol, isomericSmiles = USE_STEREOCHEMISTRY) # to canonicalize
        
        # Initialize results object
        result = RetroResult(smiles)
        
        for template in self.top_templates():
            for precursor in self.apply_one_template(mol, smiles, template, singleonly = singleonly, stop_if = stop_if):
                result.add_precursor(precursor, self.prioritizer)
                
        return result
    
    def apply_one_template(self, mol, smiles, template, singleonly = False, stop_if = False, chiral = False):
        '''
        Takes a mol object and applies a single template. Mol object should have property molAtomMapNumber:
        react_mol = clean_reactant_mapping(react_mol)
        '''
        results = []
        try:
            if template['product_smiles']:
                react_mol = Chem.MolFromSmiles(smiles + '.' + '.'.join(template['product_smiles']))
            else:
                react_mol = mol
            if chiral: 
            #RDCHIRAL still needs some testing?
                reaction = rdchiralReaction(str('(' + str(template['reaction_smarts']).replace('>>', ')>>')))
                react_mol = rdchiralReactants(smiles)
                outcomes = rdchiralRun(template['rxn'], react_mol)
            else:
                outcomes = template['rxn'].RunReactants([react_mol])
            
        except Exception as e:
            print('warning: {}'.format(e))
            print(template['reaction_smarts'])
            return []
        for j, outcome in enumerate(outcomes):
            try:
                for x in outcome:
                    x.UpdatePropertyCache()
                    Chem.SanitizeMol(x)
                    [a.SetProp('molAtomMapNumber', a.GetProp('old_molAtomMapNumber')) \
                    for a in x.GetAtoms() \
                    if 'old_molAtomMapNumber' in a.GetPropsAsDict()]
            except Exception as e:
                #print(e) # fail quietly
                continue
            smiles_list = []
            for x in outcome:
                smiles_list.extend(Chem.MolToSmiles(x, isomericSmiles = USE_STEREOCHEMISTRY).split('.'))
            if template['intra_only'] and len(smiles_list) > 1:
                continue
            
            
            precursor = RetroPrecursor(
                smiles_list = sorted(smiles_list),
                template_id = template['_id'],
                num_examples = template['count'],
                necessary_reagent = template['necessary_reagent']
                )
               
            if '.'.join(precursor.smiles_list) == smiles: continue # no transformation
            results.append(precursor)
        
        # Were we trying to stop early?
        if stop_if: 
            return False
        return results
    
    
    def dump_to_file(self, file_name):
        '''
        Write the template database to a file, of which the path in specified in the general configuration
        '''
        if not self.TEMPLATE_DB:
            MyLogger.print_and_log("No database information to output to file.", retro_transformer_loc, level = 1)
            return
            
        with open(os.path.join(gc.retro_template_data, file_name), 'wb') as file:
            pickle.dump(self.TEMPLATE_DB, file, gc.protocol)

    
    def load_from_file(self, file_name):
        '''
        Read the template database from a previously saved file, of which the path is specified in the general
        configuration
        '''
        
        if os.path.isfile(os.path.join(gc.retro_template_data, file_name)):
            with open(os.path.join(gc.retro_template_data, file_name), 'rb') as file:
                self.TEMPLATE_DB = pickle.load(file)
        else:
            MyLogger.print_and_log("No file to read data from, using online database instead.", retro_transformer_loc, level = 1)
            self.load_databases()

    def reorder(self):
        '''
        Re-orders the list of templates (self.templates) according to 
        field 'count' in descending order. This means we will apply the
        most popular templates first
        '''
        self.num_templates = len(self.templates)
        self.templates[:] = [x for x in sorted(self.templates, key = lambda z: z['count'], reverse = True)]
        self.id_to_index = {template['_id']: i for i, template in enumerate(self.templates)}

    def top_templates(self):
        '''Generator to return only top templates. 
        Assumes templates are already sorted'''
        for template in self.templates:
            if template['count'] < self.mincount: 
                break
            yield template
    def lookup_id(self, template_id):
        '''
        Find the reaction smarts for this template_id
        '''
        if template_id in self.id_to_index:
            return self.templates[self.id_to_index[template_id]]
    
    def load_databases(self):
        db_client = MongoClient(gc.MONGO['path'],gc.MONGO['id'], connect = gc.MONGO['connect'])
        self.TEMPLATE_DB = db_client[gc.RETRO_TRANSFORMS['database']][gc.RETRO_TRANSFORMS['collection']]
        
    