from __future__ import print_function
import global_config as gc
from global_config import USE_STEREOCHEMISTRY
import rdkit.Chem as Chem          
from rdkit.Chem import AllChem
import numpy as np
from functools import partial # used for passing args to multiprocessing
from utilities.i_o.logging import MyLogger
from forward_enumeration import *
from pymongo import MongoClient
from interfaces.template_transformer import TemplateTransformer
from interfaces.forward_enumerator import ForwardEnumerator
from prioritization.precursor_prioritization.heuristic_prioritizer import HeuristicPrioritizer
from prioritization.precursor_prioritization.scs_prioritizer import SCSPrioritizer
from prioritization.template_prioritization.popularity_prioritizer import PopularityPrioritizer
from prioritization.template_prioritization.relevance_prioritizer import RelevancePrioritizer
from prioritization.default_prioritizer import DefaultPrioritizer
from makeit.utilities.reactants import clean_reactant_mapping
from makeit.utilities.outcomes import summarize_reaction_outcome
from rdchiral.main import rdchiralRun
from rdchiral.initialization import rdchiralReaction, rdchiralReactants

forward_transformer_loc = 'forward_transformer'

class ForwardTransformer(TemplateTransformer, ForwardEnumerator):
    '''
    The Transformer class defines an object which can be used to perform
    one-step retrosyntheses for a given molecule.
    '''

    def __init__(self, mincount=0, TEMPLATE_DB = None, loc = False, done = None, celery = False):
        '''
        Initialize a transformer.
        TEMPLATE_DB: indicate the database you want to use (def. none)
        loc: indicate that local file data should be read instead of online data (def. false)
        '''
        
        self.done = done
        self.mincount = mincount
        self.templates = []
        self.id_to_index = {}
        self.celery = celery
        self.template_prioritizers = {}
        self.TEMPLATE_DB = TEMPLATE_DB
    
    def template_count(self):
        return len(self.templates)
    
    def get_prioritizers(self, template_prioritizer):
        if not template_prioritizer:
            MyLogger.print_and_log('Cannot run the synthetic transformer without a template prioritization method. Exiting...', forward_transformer_loc, level = 3)
        if template_prioritizer in self.template_prioritizers:
            template = self.template_prioritizers[template_prioritizer]
        else:
            if template_prioritizer == gc.popularity:
                template = PopularityPrioritizer()
            elif template_prioritizer == gc.relevance:
                template = RelevancePrioritizer()
            elif template_prioritizer == gc.natural:
                template = DefaultPrioritizer()
            else:
                template = DefaultPrioritizer()
                MyLogger.print_and_log('Prioritization method not recognized. Using natural prioritization.', forward_transformer_loc, level = 1)
                
            template.load_model()
            self.template_prioritizers[template_prioritizer] = template
    
        self.template_prioritizer = template
        
    def get_outcomes(self, smiles, mincount, template_prioritization, start_at = -1, end_at = -1, 
                     singleonly = True, stop_if = False, chiral = False):
        '''
        Each candidate in self.result.products is of type ForwardProduct
        '''
        self.get_prioritizers(template_prioritization)
        prioritized_templates = self.template_prioritizer.get_priority((self.templates, smiles))
        self.mincount = mincount
        self.start_at = start_at
        self.singleonly = singleonly
        self.stop_if = stop_if
        
        if end_at == -1 or end_at >= len(self.templates):
            self.end_at = len(self.templates)
        else:
            self.end_at = end_at
        self.chiral = chiral
         # Define mol to operate on
        mol = Chem.MolFromSmiles(smiles)
        clean_reactant_mapping(mol)
        reactants_smiles = Chem.MolToSmiles(mol)
        smiles = Chem.MolToSmiles(mol, isomericSmiles = USE_STEREOCHEMISTRY) # to canonicalize
        # Initialize results object
        if self.celery:
            result = []
        else:
            result = ForwardResult(smiles)
        
        for i in range(self.start_at,self.end_at):
            #only use templates between the specified boundaries.
            template = prioritized_templates[i]
            if template['count'] > mincount:
                products = self.apply_one_template(mol, smiles, template, singleonly=singleonly, stop_if=stop_if)
                
                if self.celery:
                    for product in products:
                        result.append({'smiles_list':product.smiles_list,
                           'smiles':product.smiles,
                           'edits':product.edits,
                           'template_ids':product.template_ids,
                           'num_examples':product.num_examples
                           })
                else:
                    result.add_products(products)
        
        return (smiles,result)
        
    def load(self, chiral = False, lowe=False, refs=False, efgs=False):
        '''
        Loads and parses the template database to a useable one
        '''
        MyLogger.print_and_log('Loading synthetic transformer, including all templates with more than {} hits'.format(self.mincount), forward_transformer_loc)
        # Save collection TEMPLATE_DB
        if not self.TEMPLATE_DB:
            self.load_databases(chiral)
        
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
            try:
                if lowe:
                    reaction_smarts_synth = '(' + reaction_smarts.split('>')[2] + ')>>(' + reaction_smarts.split('>')[0] + ')'
                else:
                    reaction_smarts_synth = '(' + reaction_smarts.replace('>>', ')>>(') + ')'
                rxn_f = AllChem.ReactionFromSmarts(reaction_smarts_synth)
                #if rxn_f.Validate() == (0, 0):
                if rxn_f.Validate()[1] == 0:
                    template['rxn_f'] = rxn_f
                else:
                    template['rxn_f'] = None
            except Exception as e:
                MyLogger.print_and_log('Couldnt load forward: {}: {}'.format(reaction_smarts_synth, e), forward_transformer_loc, level = 1)
                template['rxn_f'] = None

            if not template['rxn_f']: continue

            # Add to list
            self.templates.append(template)

        self.num_templates = len(self.templates)
        
        MyLogger.print_and_log('Synthetic transformer has been loaded - using {} templates'.format(self.num_templates), forward_transformer_loc)
        
        #multiprocessing notify done
        if self.done == None:
            pass
        else:
            self.done.value = 1
              
    
    def apply_one_template(self, mol, smiles, template, singleonly=True, stop_if=False):
        '''
        Takes a mol object and applies a single template. 
        '''

        try:
            if template['product_smiles']:
                react_mol = Chem.MolFromSmiles(smiles + '.' + '.'.join(template['product_smiles']))
            else:
                react_mol = mol
            
            if self.chiral: 
            #RDCHIRAL still needs some testing?
                reaction = rdchiralReaction(str('(' + str(template['reaction_smarts']).replace('>>', ')>>')))
                react_mol = rdchiralReactants(smiles)
                outcomes = rdchiralRun(template['rxn_f'], react_mol)
            else:
                outcomes = template['rxn_f'].RunReactants([react_mol])

        except Exception as e:
            MyLogger.print_and_log('Failed transformation for {} because of {}'.format(template['reaction_smarts'], e), forward_transformer_loc, level = 1)
            return []
        
        results = []
        if not outcomes:
            pass
        else:
            for outcome in outcomes:
                smiles_list = []
                outcome = outcome[0] # all products represented as single mol by transforms
                
                try:
                    outcome.UpdatePropertyCache()
                    Chem.SanitizeMol(outcome)
                except Exception as e:
                    MyLogger.print_and_log('Non-sensible molecule constructed by template {}'.format(template['reaction_smarts']), forward_transformer_loc, level = 1)
                    continue
                [a.SetProp(str('molAtomMapNumber'), a.GetProp(str('old_molAtomMapNumber'))) \
                    for a in outcome.GetAtoms() \
                    if str('old_molAtomMapNumber') in a.GetPropsAsDict()]
            
                # Reduce to largest (longest) product only
                candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles=True)
                smiles_list = candidate_smiles.split('.')
                if singleonly:
                    candidate_smiles = max(candidate_smiles.split('.'), key = len)
                outcome = Chem.MolFromSmiles(candidate_smiles)
                
                # Find what edits were made
                edits = summarize_reaction_outcome(react_mol, outcome)

                # Remove mapping before matching
                [x.ClearProp(str('molAtomMapNumber')) for x in outcome.GetAtoms() \
                    if x.HasProp(str('molAtomMapNumber'))] # remove atom mapping from outcome

                # Overwrite candidate_smiles without atom mapping numbers
                candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles=True)
                
                product = ForwardProduct(
                    smiles_list = sorted(smiles_list),
                    smiles = candidate_smiles,
                    template_id = str(template['_id']),
                    num_examples = template['count'],
                    edits = edits
                    )
                
                if candidate_smiles == smiles: continue # no transformation
                if stop_if:
                    if stop_if in product.smiles_list: 
                        print('Found true product - skipping remaining templates to apply')
                        return True
                else:
                    results.append(product)
            # Were we trying to stop early?
            if stop_if: 
                return False
            
        return results
    
    
    
    def dump_to_file(self, file_name):
        '''
        Write the template database to a file, of which the path in specified in the general configuration
        '''
        if not self.TEMPLATE_DB:
            MyLogger.print_and_log("No database information to output to file.", forward_transformer_loc, level = 1)
            return
    
        with open(os.path.join(gc.synth_template_data, file_name), 'wb') as file:
            pickle.dump(self.TEMPLATE_DB, file, gc.protocol)

    
    def load_from_file(self, file_name):
        '''
        Read the template database from a previously saved file, of which the path is specified in the general
        configuration
        '''

        if os.path.isfile(os.path.join(gc.synth_template_data, file_name)):
            with open(os.path.join(gc.synth_template_data, file_name), 'rb') as file:
                self.TEMPLATE_DB = pickle.load(file)
        else:
            MyLogger.print_and_log("No file to read data from, using online database instead.", forward_transformer_loc, level = 1)
            self.load_databases()

        
    def load_databases(self, chiral):
        db_client = MongoClient(gc.MONGO['path'],gc.MONGO['id'], connect = gc.MONGO['connect'])
        self.TEMPLATE_DB = db_client[gc.SYNTH_TRANSFORMS['database']][gc.SYNTH_TRANSFORMS['collection']]
    
    def top_templates(self, target):
        '''
        Generator to return only top templates. 
        First applies the template prioritization method and returns top of that list.
        '''
        prioritized_templates = self.template_prioritizer.get_priority((self.templates, target))
        counter = 0
        for template in prioritized_templates:
            # only yield template if between start and end points and if mincount criterium is fulfilled.
            # do not break on count<mincount: assumes sorted by popularity, not necessarily true
            if template['count'] < self.mincount: 
                counter += 1
            elif counter<self.start_at:
                counter += 1
            elif counter>self.end_at:
                counter += 1
            else:
                counter += 1
                yield template

            
    def lookup_id(self, template_id):
        '''
        Find the reaction smarts for this template_id
        '''
        if template_id in self.id_to_index:
            return self.templates[self.id_to_index[template_id]]     
if __name__ == '__main__':
    MyLogger.initialize_logFile()
    ft = ForwardTransformer(mincount = 100)
    ft.load()
    smiles, out = ft.get_outcomes("CCC(=O)O.CCCCN", 100, start_at=-1, end_at=-1, template_prioritization = gc.popularity)
    for outp in out.products:
        print(outp.smiles)