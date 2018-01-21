from __future__ import print_function
import makeit.global_config as gc
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from makeit.prioritization.precursors.heuristic import HeuristicPrecursorPrioritizer
from makeit.prioritization.precursors.scscore import SCScorePrecursorPrioritizer
from makeit.prioritization.templates.popularity import PopularityTemplatePrioritizer
from makeit.prioritization.templates.relevance import RelevanceTemplatePrioritizer
from makeit.prioritization.default import DefaultPrioritizer
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from pymongo import MongoClient
from makeit.utilities.io.logging import MyLogger
transformer_loc = 'template_transformer'
import cPickle as pickle
import os

class TemplateTransformer(object):
    '''
    The Transformer class defines an object which can be used to perform
    one-step retrosyntheses for a given molecule.
    '''

    def __init__(self, mincount=4, TEMPLATE_DB=None, loc=False, done=None):
        '''
        Initialize a transformer.
        TEMPLATE_DB: indicate the database you want to use (def. none)
        loc: indicate that local file data should be read instead of online data (def. false)
        '''

        # Dictionary to keep track of ID -> index in self.templates
        self.id_to_index = {}

        return

    def get_precursor_prioritizers(self, precursor_prioritizer):
        if not precursor_prioritizer:
            MyLogger.print_and_log(
                'Cannot run the Transformer without a precursor prioritization method. Exiting...', transformer_loc, level=3)
        if precursor_prioritizer in self.precursor_prioritizers:
            precursor = self.precursor_prioritizers[precursor_prioritizer]
        else:
            if precursor_prioritizer == gc.heuristic:
                precursor = HeuristicPrecursorPrioritizer()
            elif precursor_prioritizer == gc.scscore:
                precursor = SCScorePrecursorPrioritizer()
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

    def dump_to_file(self, retro, file_name, chrial=False, file_path=""):
        '''
        Write the template database to a file, of which the path in specified in the general configuration
        '''
        if not self.templates:
            self.load(chiral=chiral)

        if retro:
            if file_path=="":
                file_path = gc.retro_template_data
            if chiral:
                pickle_templates = []
                # reconstruct template list, but without chiral rxn object =>
                # can't be pickled.
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
            file = open(os.path.join(file_path, file_name), "w+")

            MyLogger.print_and_log('Wrote templates to {}'.format(
                os.path.join(file_path, file_name)), transformer_loc)

        else:
            if file_path=="":
                file_path = gc.synth_template_data
            with open(os.path.join(file_path, file_name), 'w+') as file:
                pickle.dump(self.templates, file, gc.protocol)

            MyLogger.print_and_log('Wrote templates to {}'.format(
                os.path.join(file_path, file_name)), transformer_loc)

    def load_from_file(self, retro, file_name, chiral=False, rxns=True, file_path = ""):
        '''
        Read the template database from a previously saved file, of which the path is specified in the general
        configuration
        '''
        
        MyLogger.print_and_log(
            'Loading templates from {}'.format(file_name), transformer_loc)

        if retro:
            if file_path == "":
                file_path = gc.retro_template_data
            if os.path.isfile(os.path.join(file_path, file_name)):
                with open(os.path.join(file_path, file_name), 'rb') as file:
                    if chiral and rxns:
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
                MyLogger.print_and_log(
                    "No file to read data from, using online database instead.", transformer_loc, level=1)
                self.load(chiral=chiral)

        else:
            if file_path == "":
                file_path = gc.synth_template_data
            if os.path.isfile(os.path.join(file_path, file_name)):
                with open(os.path.join(file_path, file_name), 'rb') as file:
                    self.templates = pickle.load(file)
            else:
                MyLogger.print_and_log(
                    "No file to read data from, using online database instead.", transformer_loc, level=1)
                self.load()

        self.num_templates = len(self.templates)
        MyLogger.print_and_log('Loaded templates. Using {} templates'.format(
            self.num_templates), transformer_loc)

    def get_prioritizers(self, args):
        '''
        Get the prioritization methods for the transformer (templates and/or precursors)
        '''
        raise NotImplementedError

    def load(self, chiral=False, lowe=False, refs=False, efgs=False, rxn_ex=False):
        raise NotImplementedError

    def reorder(self):
        '''Reorder self.templates in descending popularity. Also builds id_tO_index table'''
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

    def load_templates(self, retro, chiral=False, lowe=False, refs=False, rxns=True, efgs=False, rxn_ex=False):
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
                    # Force reactants and products to be one molecule (not really,
                    # but for bookkeeping)
                    if lowe:
                        reaction_smarts_one = '(' + reaction_smarts.split(
                            '>')[2] + ')>>(' + reaction_smarts.split('>')[0] + ')'
                    else:
                        reaction_smarts_one = '(' + \
                            reaction_smarts.replace('>>', ')>>(') + ')'

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
                        # if rxn_f.Validate() == (0, 0):
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

    def get_outcomes(self, smiles, mincount, prioritizers, start_at=-1, end_at=-1, singleonly=False, stop_if=False, chiral=False):
        '''
        Performs a one-step retrosynthesis given a SMILES string of a
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

    def apply_one_template(self, mol, smiles, template, singleonly=False, stop_if=False, chiral=False):
        '''
        Takes a mol object and applies a single template, returning
        a list of precursors or outcomes, depending on whether retro or 
        synthetic templates are used
        '''
        raise NotImplementedError
