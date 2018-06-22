from __future__ import print_function
import makeit.global_config as gc
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from makeit.prioritization.precursors.heuristic import HeuristicPrecursorPrioritizer
from makeit.prioritization.precursors.relevanceheuristic import RelevanceHeuristicPrecursorPrioritizer
from makeit.prioritization.precursors.mincost import MinCostPrecursorPrioritizer
from makeit.prioritization.precursors.scscore import SCScorePrecursorPrioritizer
from makeit.prioritization.templates.popularity import PopularityTemplatePrioritizer
from makeit.prioritization.templates.relevance import RelevanceTemplatePrioritizer
from makeit.prioritization.default import DefaultPrioritizer
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from pymongo import MongoClient
from makeit.utilities.io.logging import MyLogger
transformer_loc = 'template_transformer'
import makeit.utilities.io.pickle as pickle
import os

class TemplateTransformer(object):
    '''
    The Transformer class defines an object which can be used to perform
    one-step retrosyntheses for a given molecule.
    '''

    def __init__(self):
        self.id_to_index = {} # Dictionary to keep track of ID -> index in self.templates

    def get_precursor_prioritizers(self, precursor_prioritizer):
        if not precursor_prioritizer:
            MyLogger.print_and_log(
                'Cannot run the Transformer without a precursor prioritization method. Exiting...', transformer_loc, level=3)
        if precursor_prioritizer in self.precursor_prioritizers:
            precursor = self.precursor_prioritizers[precursor_prioritizer]
        else:
            if precursor_prioritizer == gc.heuristic:
                precursor = HeuristicPrecursorPrioritizer()
            elif precursor_prioritizer == gc.relevanceheuristic:
                precursor = RelevanceHeuristicPrecursorPrioritizer()
            elif precursor_prioritizer == gc.scscore:
                precursor = SCScorePrecursorPrioritizer()
            elif precursor_prioritizer == gc.mincost:
                precursor = MinCostPrecursorPrioritizer()
            elif precursor_prioritizer == gc.natural:
                precursor = DefaultPrioritizer()
            else:
                precursor = DefaultPrioritizer()
                MyLogger.print_and_log(
                    'Prioritization method not recognized. Using natural prioritization.', transformer_loc, level=1)

            precursor.load_model()
            self.precursor_prioritizers[precursor_prioritizer] = precursor

        self.precursor_prioritizer = precursor

    def get_template_prioritizers(self, template_prioritizer):
        if not template_prioritizer:
            MyLogger.print_and_log(
                'Cannot run the Transformer without a template prioritization method. Exiting...', transformer_loc, level=3)
        if template_prioritizer in self.template_prioritizers:
            template = self.template_prioritizers[template_prioritizer]
        else:
            if template_prioritizer == gc.popularity:
                template = PopularityTemplatePrioritizer()
            elif template_prioritizer == gc.relevance:
                template = RelevanceTemplatePrioritizer()
            else:
                template = PopularityTemplatePrioritizer()
                MyLogger.print_and_log('Prioritization method not recognized. Using literature popularity prioritization.', transformer_loc, level = 1)
                
            template.load_model()
            self.template_prioritizers[template_prioritizer] = template

        self.template_prioritizer = template

    def dump_to_file(self, retro, file_path, chiral=False):
        '''
        Write the template database to a file, of which the path in specified in the general configuration
        '''

        if not self.templates:
            raise ValueError('Cannot dump to file if templates have not been loaded')

        if retro and chiral:
            pickle_templates = []
            # reconstruct template list, but without chiral rxn object (can't be pickled)
            for template in self.templates:
                pickle_templates.append({
                                        'name':                 template['name'],
                                        'reaction_smarts':      template['reaction_smarts'],
                                        'incompatible_groups':  template['incompatible_groups'],
                                        'references':           template['references'],
                                        'rxn_example':          template['rxn_example'],
                                        'explicit_H':           template['explicit_H'],
                                        '_id':                  template['_id'],
                                        'product_smiles':       template['product_smiles'],
                                        'necessary_reagent':    template['necessary_reagent'],
                                        'efgs':                 template['efgs'],
                                        'intra_only':           template['intra_only'],
                                        'dimer_only':           template['dimer_only'],
                                        'chiral':               template['chiral'],
                                        'count':                template['count'],
                                        })
        else:
            pickle_templates = self.templates

        with open(file_path, 'w+') as file:
            pickle.dump(pickle_templates, file)

            MyLogger.print_and_log('Wrote templates to {}'.format(
                file_path, transformer_loc))

    def load_from_file(self, retro, file_path, chiral=False, rxns=True, refs=False, efgs=False, rxn_ex=False):
        '''
        Read the template database from a previously saved file, of which the path is specified in the general
        configuration

        retro: whether in the retrosynthetic direction
        file_path: .pickle file to read dumped templates from 
        chiral: whether to handle chirality properly (only for retro for now)
        rxns : whether or not to actually load the reaction objects (or just the info)
        '''
        
        MyLogger.print_and_log('Loading templates from {}'.format(file_path), transformer_loc)

        if os.path.isfile(file_path):
            with open(file_path, 'rb') as file:
                if retro and chiral and rxns: # cannot pickle rdchiralReactions, so need to reload from SMARTS
                    pickle_templates = pickle.load(file)
                    self.templates = []
                    for template in pickle_templates:
                        try:
                            template['rxn'] = rdchiralReaction(
                                str('(' + template['reaction_smarts'].replace('>>', ')>>(') + ')'))
                        except Exception as e:
                            template['rxn'] = None
                        self.templates.append(template)
                else:
                    self.templates = pickle.load(file)
        else:
            MyLogger.print_and_log("No file to read data from.", transformer_loc, level=1)
            raise IOError('File not found to load template_transformer from!')

        # Clear out unnecessary info
        if not refs:
            [self.templates[i].pop('refs', None) for i in range(len(self.templates))]
        elif 'refs' not in self.templates[0]:
            raise IOError('Save file does not contain references (which were requested!)')

        if not efgs:
            [self.templates[i].pop('efgs', None) for i in range(len(self.templates))]
        elif 'efgs' not in self.templates[0]:
            raise IOError('Save file does not contain efg info (which was requested!)')

        if not rxn_ex:
            [self.templates[i].pop('rxn_ex', None) for i in range(len(self.templates))]
        elif 'rxn_ex' not in self.templates[0]:
            raise IOError('Save file does not contain a reaction example (which was requested!)')


        self.num_templates = len(self.templates)
        MyLogger.print_and_log('Loaded templates. Using {} templates'.format(self.num_templates), transformer_loc)

    def get_prioritizers(self, *args, **kwargs):
        '''
        Get the prioritization methods for the transformer (templates and/or precursors)
        '''
        raise NotImplementedError

    def load(self, *args, **kwargs):
        '''
        Load and initialize templates
        '''
        raise NotImplementedError

    def reorder(self):
        '''Reorder self.templates in descending popularity. Also builds id_to_index table'''
        self.num_templates = len(self.templates)
        self.templates = sorted(self.templates, key=lambda z: z[
                                'count'], reverse=True)
        self.id_to_index = {template['_id']: i for i,
                            template in enumerate(self.templates)}
        return

    def lookup_id(self, template_id):
        '''
        Find the reaction smarts for this template_id
        '''

        if not self.id_to_index:  # need to build
            self.id_to_index = {template['_id']: i for (
                i, template) in enumerate(self.templates)}
        return self.templates[self.id_to_index[template_id]]

    def load_from_database(self, retro, chiral=False, refs=False, rxns=True, efgs=False, rxn_ex=False):
        # Save collection TEMPLATE_DB
        if not self.TEMPLATE_DB:
            self.load_databases(retro, chiral=chiral)
        self.chiral = chiral
        if self.mincount and 'count' in self.TEMPLATE_DB.find_one():
            if retro:
                filter_dict = {'count': {'$gte': min(
                    self.mincount, self.mincount_chiral)}}
            else:
                filter_dict = {'count': {'$gte': self.mincount}}
        else:
            filter_dict = {}

        # Look for all templates in collection
        to_retrieve = ['_id', 'reaction_smarts',
                       'necessary_reagent', 'count', 'intra_only', 'dimer_only']
        if refs:
            to_retrieve.append('references')
        if efgs:
            to_retrieve.append('efgs')
        if rxn_ex:
            to_retrieve.append('rxn_example')
        for document in self.TEMPLATE_DB.find(filter_dict, to_retrieve):
            # Skip if no reaction SMARTS
            if 'reaction_smarts' not in document:
                continue
            reaction_smarts = str(document['reaction_smarts'])
            if not reaction_smarts:
                continue

            if retro:
                # different thresholds for chiral and non chiral reactions
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
                'reference':            document['reference'] if 'reference' in document else '',
                'references':           document['references'] if 'references' in document else [],
                'rxn_example':          document['rxn_example'] if 'rxn_example' in document else '',
                'explicit_H':           document['explicit_H'] if 'explicit_H' in document else False,
                '_id':                  document['_id'] if '_id' in document else -1,
                'product_smiles':       document['product_smiles'] if 'product_smiles' in document else [],
                'necessary_reagent':    document['necessary_reagent'] if 'necessary_reagent' in document else '',
                'efgs':                 document['efgs'] if 'efgs' in document else None,
                'intra_only':           document['intra_only'] if 'intra_only' in document else False,
                'dimer_only':           document['dimer_only'] if 'dimer_only' in document else False,
            }
            if retro:
                template['chiral'] = chiral_rxn

            # Frequency/popularity score
            if 'count' in document:
                template['count'] = document['count']
            elif 'popularity' in document:
                template['count'] = document['popularity']
            else:
                template['count'] = 1

            # Define reaction in RDKit and validate
            if rxns:
                try:
                    # Force reactants and products to be one pseudo-molecule (bookkeeping)
                    reaction_smarts_one = '(' + reaction_smarts.replace('>>', ')>>(') + ')'

                    if retro:
                        if chiral:
                            rxn = rdchiralReaction(str(reaction_smarts_one))
                            template['rxn'] = rxn
                        else:
                            rxn = AllChem.ReactionFromSmarts(
                                str(reaction_smarts_one))
                            if rxn.Validate()[1] == 0:
                                template['rxn'] = rxn
                            else:
                                template['rxn'] = None
                    else:
                        rxn_f = AllChem.ReactionFromSmarts(reaction_smarts_one)
                        if rxn_f.Validate()[1] == 0:
                            template['rxn_f'] = rxn_f
                        else:
                            template['rxn_f'] = None

                except Exception as e:
                    if gc.DEBUG:
                        MyLogger.print_and_log('Couldnt load : {}: {}'.format(
                            reaction_smarts_one, e), transformer_loc, level=1)
                    template['rxn'] = None
                    template['rxn_f'] = None

            # Add to list
            self.templates.append(template)

        self.reorder()

    def get_outcomes(self, *args, **kwargs):
        '''
        Performs a one-step transformation given a SMILES string of a
        target molecule by applying each transformation template
        sequentially.
        '''
        raise NotImplementedError

    def load_databases(self, retro, chiral=False):

        db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                                'id'], connect=gc.MONGO['connect'])

        if retro:
            if self.chiral:
                self.TEMPLATE_DB = db_client[gc.RETRO_TRANSFORMS_CHIRAL[
                    'database']][gc.RETRO_TRANSFORMS_CHIRAL['collection']]
                MyLogger.print_and_log("Using {} as template database.".format(
                    gc.RETRO_TRANSFORMS_CHIRAL['collection']), retro_transformer_loc)
            else:
                self.TEMPLATE_DB = db_client[gc.RETRO_TRANSFORMS[
                    'database']][gc.RETRO_TRANSFORMS['collection']]
                MyLogger.print_and_log("Using {} as template database.".format(
                    gc.RETRO_TRANSFORMS['collection']), retro_transformer_loc)
        else:
            self.TEMPLATE_DB = db_client[gc.SYNTH_TRANSFORMS[
                'database']][gc.SYNTH_TRANSFORMS['collection']]

    def apply_one_template(self, *args, **kwargs):
        '''
        Takes a mol object and applies a single template, returning
        a list of precursors or outcomes, depending on whether retro or 
        synthetic templates are used
        '''
        raise NotImplementedError
