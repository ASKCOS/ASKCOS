from __future__ import print_function

import makeit.global_config as gc
import os
import cPickle as pickle
from pymongo import MongoClient

USE_STEREOCHEMISTRY = True
import rdkit.Chem as Chem          
from rdkit.Chem import AllChem
import numpy as np
from functools import partial # used for passing args to multiprocessing
from makeit.utilities.i_o.logging import MyLogger
from makeit.utilities.reactants import clean_reactant_mapping
from makeit.retro_synthetic.results import RetroResult, RetroPrecursor
from makeit.interfaces.template_transformer import TemplateTransformer
from makeit.prioritization.precursors.heuristic import HeuristicPrecursorPrioritizer
from makeit.prioritization.precursors.scscore import SCScorePrecursorPrioritizer
from makeit.prioritization.templates.popularity import PopularityTemplatePrioritizer
from makeit.prioritization.templates.relevance import RelevanceTemplatePrioritizer
from makeit.prioritization.default import DefaultPrioritizer
from rdchiral.main import rdchiralRun
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
retro_transformer_loc = 'retro_transformer'

class RetroTransformer(TemplateTransformer):
    '''
    The Transformer class defines an object which can be used to perform
    one-step retrosyntheses for a given molecule.
    '''

    def __init__(self, celery = False, mincount=0, mincount_chiral=-1, TEMPLATE_DB = None, loc = False, done = None):
        '''
        Initialize a transformer.
        TEMPLATE_DB: indicate the database you want to use (def. none)
        loc: indicate that local file data should be read instead of online data (def. false)
        '''
        self.done = done
        self.mincount = mincount
        if mincount_chiral == -1:
            self.mincount_chiral = mincount
        else:
            self.mincount_chiral = mincount_chiral
        self.templates = []
        self.celery = celery
        #if not self.celery and not TEMPLATE_DB:
        #    MyLogger.print_and_log('Predefined template database is required for the retro transformer. Exiting...', retro_transformer_loc, level=3)
        
        self.TEMPLATE_DB = TEMPLATE_DB
        self.precursor_prioritizers = {}
        self.template_prioritizers = {}
        self.precursor_prioritizer = None
        self.template_prioritizer = None
        
        super(RetroTransformer, self).__init__()
        
    def get_prioritizers(self, precursor_prioritizer = None, template_prioritizer = None):
        
        if not precursor_prioritizer:
            MyLogger.print_and_log('Cannot run the Retro-Transformer without a precursor prioritization method. Exiting...', retro_transformer_loc, level = 3)
        if precursor_prioritizer in self.precursor_prioritizers:
            precursor = self.precursor_prioritizers[precursor_prioritizer]
        else:
            if precursor_prioritizer == gc.heuristic:
                precursor = HeuristicPrioritizer()
            elif precursor_prioritizer == gc.scscore:
                precursor = SCScorePrioritizer()
            elif precursor_prioritizer == gc.natural:
                precursor = DefaultPrioritizer()
            else:
                precursor = DefaultPrioritizer()
                MyLogger.print_and_log('Prioritization method not recognized. Using natural prioritization.', retro_transformer_loc, level = 1)
                
            precursor.load_model()
            self.precursor_prioritizers[precursor_prioritizer] = precursor
        
        if not template_prioritizer:
            MyLogger.print_and_log('Cannot run the Retro-Transformer without a template prioritization method. Exiting...', retro_transformer_loc, level = 3)
        if template_prioritizer in self.template_prioritizers:
            template = self.template_prioritizers[template_prioritizer]
        else:
            if template_prioritizer == gc.popularity:
                template = PopularityPrioritizer()
            elif template_prioritizer == gc.relevance:
                template = RelevancePrioritizer()
            else:
                template = PopularityPrioritizer()
                MyLogger.print_and_log('Prioritization method not recognized. Using literature popularity prioritization.', retro_transformer_loc, level = 1)
                
            template.load_model()
            self.template_prioritizers[template_prioritizer] = template
    
        self.precursor_prioritizer = precursor
        self.template_prioritizer = template
        
    def load(self, TEMPLATE_DB=None, chiral=False, lowe=False, refs=False, rxns=True, efgs=False, rxn_ex=False):
        '''
        Loads and parses the template database to a useable one
        '''
       
        MyLogger.print_and_log('Loading retro-synthetic transformer, including all templates with more than {} hits ({} for chiral reactions)'.format(self.mincount, self.mincount_chiral), retro_transformer_loc)
        self.chiral = chiral

        # Override previous TEMPLATE_DB?
        if TEMPLATE_DB is not None:
            self.TEMPLATE_DB = TEMPLATE_DB
        
        # Save collection TEMPLATE_DB
        if not self.TEMPLATE_DB:
            self.load_databases()
        
        if self.mincount and 'count' in self.TEMPLATE_DB.find_one(): 
            filter_dict = {'count': { '$gte': min(self.mincount,self.mincount_chiral)}}
        else: 
            filter_dict = {}
       
        # Look for all templates in collection
        to_retrieve = ['_id', 'reaction_smarts', 'necessary_reagent', 'count', 'intra_only','dimer_only']
        if refs:
            to_retrieve.append('references')
        if efgs:
            to_retrieve.append('efgs')
        if rxn_ex:
            to_retrieve.append('rxn_example')
        for document in self.TEMPLATE_DB.find(filter_dict, to_retrieve):
            # Skip if no reaction SMARTS
            if 'reaction_smarts' not in document: continue
            reaction_smarts = str(document['reaction_smarts'])
            if not reaction_smarts: continue
            
            #different thresholds for chiral and non chiral reactions
            chiral_rxn = False
            for c in reaction_smarts:
                if c in ('@', '/', '\\'):
                    chiral_rxn = True 
                    break

            if chiral_rxn and document['count'] < self.mincount_chiral:
                continue
            if not chiral_rxn and document['count'] < self.mincount:
                continue
            
            
            # Define dictionary
            template = {
                'name':                 document['name'] if 'name' in document else '',
                'reaction_smarts':      reaction_smarts,
                'incompatible_groups':  document['incompatible_groups'] if 'incompatible_groups' in document else [],
                'references':           document['references'] if 'references' in document else [],
                'rxn_example':          document['rxn_example'] if 'rxn_example' in document else '',
                'explicit_H':           document['explicit_H'] if 'explicit_H' in document else False,
                '_id':                  document['_id'] if '_id' in document else -1,
                'product_smiles':       document['product_smiles'] if 'product_smiles' in document else [], 
                'necessary_reagent':    document['necessary_reagent'] if 'necessary_reagent' in document else '',       
                'efgs':                 document['efgs'] if 'efgs' in document else None,
                'intra_only':           document['intra_only'] if 'intra_only' in document else False,
                'dimer_only':           document['dimer_only'] if 'dimer_only' in document else False,        
                'chiral':               chiral_rxn
            }

            # Frequency/popularity score
            if 'count' in document: 
                template['count'] = document['count']
            elif 'popularity' in document:
                template['count'] = document['popularity']
            else:
                template['count'] = 1
            
            if rxns:  # Define reaction in RDKit and validate
                try:
                    # Force reactants and products to be one molecule (not really, but for bookkeeping)
                    reaction_smarts_retro = '(' + reaction_smarts.replace('>>', ')>>(') + ')'
                
                    if chiral:
                        rxn = rdchiralReaction(str(reaction_smarts_retro))
                        template['rxn'] = rxn
                    else:
                        rxn = AllChem.ReactionFromSmarts(str(reaction_smarts_retro))
                        if rxn.Validate()[1] == 0: 
                            template['rxn'] = rxn
                        else:
                            template['rxn'] = None
                
                except Exception as e:
                    if gc.DEBUG:
                        MyLogger.print_and_log('Couldnt load retro: {}: {}'.format(reaction_smarts_retro, e),retro_transformer_loc, level=1)
                    template['rxn'] = None
            
            # Add to list
            self.templates.append(template)

        self.num_templates = len(self.templates)
        
        self.reorder()        
        MyLogger.print_and_log('Retro-synthetic transformer has been loaded - using {} templates.'.format(self.num_templates), retro_transformer_loc)
        
    def get_outcomes(self, smiles, mincount, prioritizers, start_at = -1, end_at = -1,
                     singleonly = False, stop_if = False):
        '''
        Performs a one-step retrosynthesis given a SMILES string of a
        target molecule by applying each transformation template
        sequentially.
        '''
        (precursor_prioritizer, template_prioritizer) = prioritizers
        #Check modules:
        if not (template_prioritizer and precursor_prioritizer):
            MyLogger.print_and_log('Template prioritizer and/or precursor prioritizer are missing. Exiting...', retro_transformer_loc, level = 3)
        self.mincount = mincount
        self.get_prioritizers(precursor_prioritizer, template_prioritizer)
        # Define mol to operate on
        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True) # to canonicalize
        if self.chiral:
            mol = rdchiralReactants(smiles)
        # Initialize results object
        result = RetroResult(smiles)
        
        for template in self.top_templates(smiles):
            for precursor in self.apply_one_template(mol, smiles, template, singleonly = singleonly, stop_if = stop_if):
                result.add_precursor(precursor, self.precursor_prioritizer)
        
        return result
    
    def apply_one_template(self, react_mol, smiles, template, singleonly = False, stop_if = False):
        '''
        Takes a mol object and applies a single template. Mol object should have property molAtomMapNumber:
        react_mol = clean_reactant_mapping(react_mol)
        '''
        results = []
        try:
            if self.chiral:
                outcomes = rdchiralRun(template['rxn'], react_mol)
            else:
                outcomes = template['rxn'].RunReactants([react_mol])
            
        except Exception as e:
            #MyLogger.print_and_log('{}'.format(e), retro_transformer_loc, level = 1)
            #print(template['reaction_smarts'])
            return []
        for j, outcome in enumerate(outcomes):
            smiles_list = []
            #Output of rdchiral is (a list of) smiles of the products.
            if self.chiral:
                smiles_list = outcome.split('.')
            #Output of the standard reactor in rdkit is an rdkit molecule object.
            else:
                try:
                    for x in outcome:
                        x.UpdatePropertyCache()
                        Chem.SanitizeMol(x)
                        [a.SetProp('molAtomMapNumber', a.GetProp('old_molAtomMapNumber')) \
                        for a in x.GetAtoms() \
                        if 'old_molAtomMapNumber' in a.GetPropsAsDict()]
                        smiles_list.extend(Chem.MolToSmiles(x, isomericSmiles = USE_STEREOCHEMISTRY).split('.'))
                except Exception as e:
                    #print(e) # fail quietly
                    continue
                
            if template['intra_only'] and len(smiles_list) > 1:
                #Disallowed intermolecular reaction
                continue
            if template['dimer_only'] and (len(set(smiles_list)) != 1 or len(smiles_list) != 2):
                #Not a dimer
                continue
            if '.'.join(smiles_list) == smiles:
                # no transformation
                continue 
            
            precursor = RetroPrecursor(
                smiles_list = sorted(smiles_list),
                template_id = str(template['_id']),
                template_score = template['score'],
                num_examples = template['count'],
                necessary_reagent = template['necessary_reagent']
                )
               
            results.append(precursor)
        
        # Were we trying to stop early?
        if stop_if: 
            return False
        return results
    
    
    def dump_to_file(self, file_name, chiral=True):
        '''
        Write the template database to a file, of which the path in specified in the general configuration
        '''
        if not self.templates:
            raise ValueError('RetroTransformer cannot dump_to_file before tempaltes initialized')
        
        if chiral: # cannot pickle pre-loaded rxn
            pickle_templates = []
            for template in self.templates:
                del template['rxn']
                pickle_templates.append(template)
        else:
            pickle_templates = self.templates
        file = open(os.path.join(gc.retro_template_data, file_name), "w+")
        
        pickle.dump(pickle_templates, file, gc.protocol)
        
        MyLogger.print_and_log('Wrote templates to {}'.format(os.path.join(gc.retro_template_data, file_name)), retro_transformer_loc)

    
    def load_from_file(self, file_name, chiral=True, rxns=True):
        '''
        Read the template database from a previously saved file, of which the path is specified in the general
        configuration
        '''
        MyLogger.print_and_log('Loading templates from {}'.format(os.path.join(gc.retro_template_data, file_name)), retro_transformer_loc)
        if os.path.isfile(os.path.join(gc.retro_template_data, file_name)):
            with open(os.path.join(gc.retro_template_data, file_name), 'rb') as file:
                if chiral and rxns:
                    pickle_templates = pickle.load(file)
                    self.templates = []
                    for template in pickle_templates:
                        try:
                            template['rxn'] = rdchiralReaction(str('(' + template['reaction_smarts'].replace('>>', ')>>(') + ')'))
                        except Exception as e:
                            template['rxn'] = None
                        self.templates.append(template)
                else:
                    self.templates = pickle.load(file)
        else:
            MyLogger.print_and_log("No file to read data from, using online database instead.", retro_transformer_loc, level = 1)
            self.load(chiral=chiral)
        self.num_templates = len(self.templates)
        MyLogger.print_and_log('Loaded templates. Using {} templates'.format(self.num_templates), retro_transformer_loc)

    def top_templates(self, target):
        '''
        Generator to return only top templates. 
        First applies the template prioritization method and returns top of that list.
        '''
        prioritized_templates = self.template_prioritizer.get_priority((self.templates, target))
        for template in prioritized_templates:
            if template['count'] < self.mincount: 
                pass
            else:
                yield template
    
    def load_databases(self):
        
        db_client = MongoClient(gc.MONGO['path'],gc.MONGO['id'], connect = gc.MONGO['connect'])
        if self.chiral:
            self.TEMPLATE_DB = db_client[gc.RETRO_TRANSFORMS_CHIRAL['database']][gc.RETRO_TRANSFORMS_CHIRAL['collection']]
            MyLogger.print_and_log("Using {} as template database.".format(gc.RETRO_TRANSFORMS_CHIRAL['collection']), retro_transformer_loc)
        else:
            self.TEMPLATE_DB = db_client[gc.RETRO_TRANSFORMS['database']][gc.RETRO_TRANSFORMS['collection']]
            MyLogger.print_and_log("Using {} as template database.".format(gc.RETRO_TRANSFORMS['collection']), retro_transformer_loc)
            
if __name__ == '__main__':
    MyLogger.initialize_logFile()
    db_client = MongoClient(gc.MONGO['path'],gc.MONGO['id'], connect = gc.MONGO['connect'])
    TEMPLATE_DB = db_client[gc.RETRO_TRANSFORMS_CHIRAL['database']][gc.RETRO_TRANSFORMS_CHIRAL['collection']]
    t = RetroTransformer(mincount = 100, mincount_chiral = 50, TEMPLATE_DB=TEMPLATE_DB)
    t.load(chiral = True)
    print(t.get_outcomes('C1C(=O)OCC12CC(C)CC2', 100, (gc.heuristic, gc.popularity)).precursors)
    
