from __future__ import print_function

import makeit.global_config as gc
import os
import six; from six.moves import cPickle as pickle
from pymongo import MongoClient

USE_STEREOCHEMISTRY = True
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from functools import partial  # used for passing args to multiprocessing
from makeit.utilities.io.logging import MyLogger
from makeit.utilities.reactants import clean_reactant_mapping
from makeit.retrosynthetic.results import RetroResult, RetroPrecursor
from makeit.interfaces.template_transformer import TemplateTransformer
from makeit.prioritization.precursors.heuristic import HeuristicPrecursorPrioritizer
from makeit.prioritization.precursors.relevanceheuristic import RelevanceHeuristicPrecursorPrioritizer
from makeit.prioritization.precursors.mincost import MinCostPrecursorPrioritizer
from makeit.prioritization.precursors.scscore import SCScorePrecursorPrioritizer
from makeit.prioritization.templates.popularity import PopularityTemplatePrioritizer
from makeit.prioritization.templates.relevance import RelevanceTemplatePrioritizer
from makeit.prioritization.default import DefaultPrioritizer
from makeit.synthetic.evaluation.fast_filter import FastFilterScorer
from rdchiral.main import rdchiralRun
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
retro_transformer_loc = 'retro_transformer'


class RetroTransformer(TemplateTransformer):
    """
    The Transformer class defines an object which can be used to perform
    one-step retrosyntheses for a given molecule.
    """

    def __init__(self, celery=False, mincount=gc.RETRO_TRANSFORMS_CHIRAL['mincount'], 
        mincount_chiral=gc.RETRO_TRANSFORMS_CHIRAL['mincount_chiral'], 
        TEMPLATE_DB=None):
        """Initialize
        
        Keyword Arguments:
            celery {bool} -- Whether or not Celery is being used (default: {False})
            mincount {number} -- Minimum number of precedents for an achiral template
                for inclusion in the template library. Only used when retrotransformers
                need to be initialized (default: {25})
            mincount_chiral {number} -- Minimum number of precedents for a chiral template
                for inclusion in the template library. Only used when retrotransformers
                need to be initialized. Chiral templates are necessarily more specific,
                so we generally use a lower threshold than achiral templates (default: {10})
            TEMPLATE_DB {None or MongoDB} -- Database to load templates from (default: {None})
        """
        
        self.mincount = mincount
        if mincount_chiral == -1:
            self.mincount_chiral = mincount
        else:
            self.mincount_chiral = mincount_chiral
        self.templates = []
        self.celery = celery
        self.TEMPLATE_DB = TEMPLATE_DB
        self.precursor_prioritizers = {}
        self.template_prioritizers = {}
        self.precursor_prioritizer = None
        self.template_prioritizer = None
        self.fast_filter = None
        if self.celery:
            # Pre-load fast filter
            self.load_fast_filter()

        super(RetroTransformer, self).__init__()

    def load(self, chiral=True, refs=False, rxns=True, efgs=False, rxn_ex=False):
        """Load templates to finish initializing the transformer
        
        Keyword Arguments:
            TEMPLATE_DB {None or MongoDB} -- MongoDB to load from (default: {None})
            chiral {bool} -- Whether to pay close attention to chirality (default: {False})
            refs {bool} -- Whether to also save references (Reaxys instance IDs)
                when loading templates (default: {False})
            rxns {bool} -- Whether to actually load reaction SMARTS into RDKit
                reaction objects (default: {True})
            efgs {bool} -- Whether to load statistics about DFG popularity [old]
                (default: {False})
            rxn_ex {bool} -- Whether to also save a reaction example with
                each template as it is loaded (default: {False})
        """


        MyLogger.print_and_log('Loading retro-synthetic transformer, including all templates with more than {} hits ({} for chiral reactions)'.format(
            self.mincount, self.mincount_chiral), retro_transformer_loc)

        if chiral:
            from makeit.utilities.io.files import get_retrotransformer_chiral_path
            file_path = get_retrotransformer_chiral_path(
                gc.RETRO_TRANSFORMS_CHIRAL['database'],
                gc.RETRO_TRANSFORMS_CHIRAL['collection'],
                gc.RETRO_TRANSFORMS_CHIRAL['mincount'],
                gc.RETRO_TRANSFORMS_CHIRAL['mincount_chiral'],
            )
        else:
            from makeit.utilities.io.files import get_retrotransformer_achiral_path
            file_path = get_retrotransformer_achiral_path(
                gc.RETRO_TRANSFORMS['database'],
                gc.RETRO_TRANSFORMS['collection'],
                gc.RETRO_TRANSFORMS['mincount'],
            )

        try:
            self.load_from_file(True, file_path, chiral=chiral, rxns=rxns, refs=refs, efgs=efgs, rxn_ex=rxn_ex)
        except IOError:
            self.load_from_database(True, chiral=chiral, rxns=rxns, refs=refs, efgs=efgs, rxn_ex=rxn_ex)
        finally:
            self.reorder()

        MyLogger.print_and_log('Retrosynthetic transformer has been loaded - using {} templates.'.format(
            self.num_templates), retro_transformer_loc)

    def load_fast_filter(self):
        # NOTE: Keras backend must be Theano for fast filter to work
        self.fast_filter = FastFilterScorer()
        # self.fast_filter.set_keras_backend('theano')
        self.fast_filter.load(model_path=gc.FAST_FILTER_MODEL['trained_model_path'])

    def get_outcomes(self, smiles, mincount, prioritizers, **kwargs):
        """Performs a one-step retrosynthesis given a SMILES string of a
        target molecule by applying each transformation template
        sequentially.
        
        Arguments:
            smiles {string} -- product SMILES string to find precursors for
            mincount {int} -- Minimum template popularity
            prioritizers {2-tuple of (string, string)} -- tuple defining the
                precursor_prioritizer and template_prioritizer to use for 
                expansion, each as a string
            **kwargs -- Additional kwargs to pass through to prioritizers or to
                handle deprecated options            
        
        Returns:
             RetroResult -- special object for a retrosynthetic expansion result,
                defined by ./results.py
        """
        apply_fast_filter = kwargs.pop('apply_fast_filter', True)
        filter_threshold = kwargs.pop('filter_threshold', 0.75)
        if (apply_fast_filter and not self.fast_filter):
            self.load_fast_filter()

        (precursor_prioritizer, template_prioritizer) = prioritizers
        # Check modules:
        if not (template_prioritizer and precursor_prioritizer):
            MyLogger.print_and_log(
                'Template prioritizer and/or precursor prioritizer are missing. Exiting...', retro_transformer_loc, level=3)
        self.mincount = mincount
        self.get_precursor_prioritizers(precursor_prioritizer)
        self.get_template_prioritizers(template_prioritizer)

        # Define mol to operate on
        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)  # to canonicalize
        if self.chiral:
            mol = rdchiralReactants(smiles)

        # Initialize results object
        result = RetroResult(smiles)

        for template in self.top_templates(smiles, **kwargs):
            for precursor in self.apply_one_template(mol, smiles, template):
                # Should we add this to the results?
                if apply_fast_filter:
                    reactant_smiles = '.'.join(precursor.smiles_list)
                    filter_flag = self.fast_filter.filter_with_threshold(reactant_smiles, smiles, filter_threshold)
                    if filter_flag:
                        result.add_precursor(precursor, self.precursor_prioritizer, **kwargs)
                else:
                    result.add_precursor(precursor, self.precursor_prioritizer, **kwargs)
                
        return result

    def apply_one_template(self, react_mol, smiles, template, **kwargs):
        """Takes a mol object and applies a single template
                
        Arguments:
            react_mol {rdchiralReactants} -- Initialized reactant object using
                RDChiral helper package; is the target compound to find 
                precursors for
            smiles {string} -- Product SMILES (no atom mapping)
            template {dict} -- Template to be applied, containing an initialized
                rdchiralReaction object as its 'rxn' field
            **kwargs -- Additional kwargs to accept deprecated options
        
        Returns:
            list -- list of RetroPrecursor objects resulting from applying
                this one template
        """
        results = []
        try:
            if self.chiral:
                outcomes = rdchiralRun(template['rxn'], react_mol)
            else:
                outcomes = template['rxn'].RunReactants([react_mol])
        except Exception as e:
            # print(template['reaction_smarts'])
            return []
        for j, outcome in enumerate(outcomes):
            smiles_list = []
            # Output of rdchiral is (a list of) smiles of the products.
            if self.chiral:
                smiles_list = outcome.split('.')
            # Output of the standard reactor in rdkit is an rdkit molecule object.
            else:
                try:
                    for x in outcome:
                        x.UpdatePropertyCache()
                        Chem.SanitizeMol(x)
                        [a.SetProp('molAtomMapNumber', a.GetProp('old_molAtomMapNumber'))
                         for a in x.GetAtoms()
                         if 'old_molAtomMapNumber' in a.GetPropsAsDict()]
                        smiles_list.extend(Chem.MolToSmiles(
                            x, isomericSmiles=USE_STEREOCHEMISTRY).split('.'))
                except Exception as e:
                    # print(e) # fail quietly
                    continue

            if template['intra_only'] and len(smiles_list) > 1:
                # Disallowed intermolecular reaction
                continue
            if template['dimer_only'] and (len(set(smiles_list)) != 1 or len(smiles_list) != 2):
                # Not a dimer
                continue
            if '.'.join(smiles_list) == smiles:
                # no transformation
                continue

            precursor = RetroPrecursor(
                smiles_list=sorted(smiles_list),
                template_id=str(template['_id']),
                template_score=template['score'],
                num_examples=template['count'],
                necessary_reagent=template['necessary_reagent']
            )

            results.append(precursor)

        return results

    def top_templates(self, target, **kwargs):
        """Generator to return only top templates. 
                
        Arguments:
            target {string} -- SMILES string of target product
            **kwargs -- additional options to pass template_prioritizer
        
        Yields:
            dict -- single templates in order of decreasing priority
        """
        prioritized_templates = self.template_prioritizer.get_priority((self.templates, target), **kwargs)
        for template in prioritized_templates:
            if template['count'] < self.mincount:
                pass
            else:
                yield template

if __name__ == '__main__':
    MyLogger.initialize_logFile()
    t = RetroTransformer()
    t.load(chiral=True)

    print(t.get_outcomes('C1C(=O)OCC12CC(C)CC2', 100,
                         (gc.relevanceheuristic, gc.relevance), max_cum_prob = 0.5).precursors)