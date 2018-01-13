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
from makeit.utilities.io.logging import MyLogger
from makeit.utilities.reactants import clean_reactant_mapping
from makeit.retrosynthetic.results import RetroResult, RetroPrecursor
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
   
    def load(self, TEMPLATE_DB=None, chiral=False, lowe=False, refs=False, rxns=True, efgs=False, rxn_ex=False):
        '''
        Loads and parses the template database to a useable one
        '''
  
        MyLogger.print_and_log('Loading retro-synthetic transformer, including all templates with more than {} hits ({} for chiral reactions)'.format(self.mincount, self.mincount_chiral), retro_transformer_loc)
   
        self.load_templates(True, chiral=chiral, lowe=lowe, refs=refs, efgs=efgs,rxn_ex = rxn_ex)
        
        MyLogger.print_and_log('Retro-synthetic transformer has been loaded - using {} templates.'.format(self.num_templates), retro_transformer_loc)
        
    def get_outcomes(self, smiles, mincount, prioritizers, start_at = -1, end_at = -1,
                     singleonly = False, stop_if = False, template_count = 10000):
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
        self.template_count = template_count
        self.get_precursor_prioritizers(precursor_prioritizer)
        self.get_template_prioritizers(template_prioritizer)
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
    
    def top_templates(self, target):
        '''
        Generator to return only top templates. 
        First applies the template prioritization method and returns top of that list.
        '''
        prioritized_templates = self.template_prioritizer.get_priority((self.templates, target), template_count = self.template_count)
        for template in prioritized_templates:
            if template['count'] < self.mincount: 
                pass
            else:
                yield template
        
if __name__ == '__main__':
    MyLogger.initialize_logFile()
    db_client = MongoClient(gc.MONGO['path'],gc.MONGO['id'], connect = gc.MONGO['connect'])
    TEMPLATE_DB = db_client[gc.RETRO_TRANSFORMS_CHIRAL['database']][gc.RETRO_TRANSFORMS_CHIRAL['collection']]
    t = RetroTransformer(mincount = 100, mincount_chiral = 50, TEMPLATE_DB=TEMPLATE_DB)
    t.load(chiral = True)
    print(t.get_outcomes('C1C(=O)OCC12CC(C)CC2', 100, (gc.heuristic, gc.popularity)).precursors)
    
