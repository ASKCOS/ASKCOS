import gzip
import json
import makeit.global_config as gc
import pymongo
from rdkit.Chem import AllChem
from rdchiral.initialization import rdchiralReaction
from pymongo import MongoClient
from bson.objectid import ObjectId
from makeit.utilities.io.logger import MyLogger
transformer_loc = 'template_transformer'
import makeit.utilities.io.pickle as pickle
import os, sys

class TemplateTransformer(object):
    """One-step retrosynthesis transformer.

    The TemplateTransformer class defines an object which can be used to perform
    one-step retrosyntheses for a given molecule.

    Attributes:
        id_to_index (dict): Maps template ID to index in ``self.templates``.
        precursor_prioritizers ():
        precursor_prioritizer ():
        template_prioritizers ():
        template_prioritizer ():
        templates ():
        num_templates (int): Number of templates loaded by the transformer.
        chiral (bool): Whether to properly handle chirality.
        mincount (int): Minimum template popularity.
        mincount_chiral (int): Minimum template popularity for chiral templates.
        TEMPLATE_DB ():
    """

    def __init__(self, load_all=gc.PRELOAD_TEMPLATES, use_db=True, TEMPLATE_DB=None):
        """Initializes TemplateTransformer.

        Args:
            load_all (bool, optional): Whether to load all of the templates into
                memory. (default: {gc.PRELOAD_TEMPLATES})
            use_db (bool, optional): Whether to use the database to look up
                templates. (default: {True})
        """
        self.templates = []
        self.load_all = load_all
        self.use_db = use_db
        self.TEMPLATE_DB = TEMPLATE_DB
        self.id_to_index = {} # Dictionary to keep track of ID -> index in self.templates

    def doc_to_template(self, document, retro=True):
        """Returns a template given a document from the database or file.

        Args:
            document (dict): Document of template from database or file.

        Returns:
            dict: Retrosynthetic template.
        """
        if 'reaction_smarts' not in document:
            return
        reaction_smarts = str(document['reaction_smarts'])
        if not reaction_smarts:
            return

        if not retro:
            document['rxn_f'] = AllChem.ReactionFromSmarts(reaction_smarts)
            return document

        # different thresholds for chiral and non chiral reactions
        chiral_rxn = False
        for c in reaction_smarts:
            if c in ('@', '/', '\\'):
                chiral_rxn = True
                break

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
        template['chiral'] = chiral_rxn

        # Frequency/popularity score
        template['count'] = document.get('count', 1)

        # Define reaction in RDKit and validate
        try:
            # Force reactants and products to be one pseudo-molecule (bookkeeping)
            reaction_smarts_one = '(' + reaction_smarts.replace('>>', ')>>(') + ')'

            rxn = rdchiralReaction(str(reaction_smarts_one))
            template['rxn'] = rxn

        except Exception as e:
            if gc.DEBUG:
                MyLogger.print_and_log('Couldnt load : {}: {}'.format(
                    reaction_smarts_one, e), transformer_loc, level=1)
            template['rxn'] = None
        return template

    def dump_to_file(self, retro, file_path, chiral=False):
        """Write the template database to a file.

        Args:
            retro (bool): Whether in the retrosynthetic direction.
            file_path (str): Specifies where to save the database.
            chiral (bool, optional): Whether to care about chirality.
                (default: {False})
        """

        if self.use_db:
            MyLogger.print_and_log('Cannot write templates when using db', transformer_loc)
            return
        if not self.templates:
            raise ValueError('Cannot dump to file if templates have not been loaded')
        
        templates = []

        if retro and chiral:
            # reconstruct template list, but without chiral rxn object
            for template in self.templates:
                templates.append({
                    'name': template['name'],
                    'reaction_smarts': template['reaction_smarts'],
                    'incompatible_groups': template['incompatible_groups'],
                    'references': template['references'],
                    'rxn_example': template['rxn_example'],
                    'explicit_H': template['explicit_H'],
                    '_id': template['_id'],
                    'product_smiles': template['product_smiles'],
                    'necessary_reagent': template['necessary_reagent'],
                    'efgs': template['efgs'],
                    'intra_only': template['intra_only'],
                    'dimer_only': template['dimer_only'],
                    'chiral': template['chiral'],
                    'count': template['count'],
                })
        else:
            templates = self.templates

        if file_path[-2:] != 'gz':
            file_path += '.gz'

        with gzip.open(file_path, 'wb') as f:
            json.dump(templates, f)

        MyLogger.print_and_log('Wrote templates to {}'.format(file_path), transformer_loc)

    def load_from_file(self, file_path, template_set=None, retro=True):
        """Read the template database from a previously saved file.

        Args:
            file_path (str): gzipped json file to read dumped templates from.
            template_set (str): optional name of template set to load, otherwisse load templates from all template sets in file
            retro (bool): whether or not templates being loaded represent retrsynthetic templates
        """

        MyLogger.print_and_log('Loading templates from {}'.format(file_path), transformer_loc)

        if os.path.isfile(file_path):
            with gzip.open(file_path, 'rb') as f:
                self.templates = json.loads(f.read().decode('utf-8'))
        else:
            MyLogger.print_and_log("No file to read data from.", transformer_loc, level=1)
            raise IOError('File not found to load template_transformer from!')

        if template_set is not None and template_set != 'all':
            self.templates = list(
                filter(
                    lambda x: x.get('template_set') == template_set, 
                    self.templates
                )
            )
        
        for n, template in enumerate(self.templates):
            if self.load_all:
                template = self.doc_to_template(template, retro=retro)
                self.templates[n] = template
            if template.get('_id') is None:
                template['_id'] = n
            self.id_to_index[template.get('_id')] = n

        self.num_templates = len(self.templates)
        MyLogger.print_and_log('Loaded templates. Using {} templates'.format(self.num_templates), transformer_loc)

    def get_prioritizers(self, *args, **kwargs):
        """Get the prioritization methods for the transformer."""
        raise NotImplementedError

    def load(self, *args, **kwargs):
        """Load and initialize templates."""
        raise NotImplementedError

    def lookup_id(self, template_id):
        """Find the reaction SMARTS for this template_id.

        Args:
            template_id (str, bytes, or ObjectId): ID of requested template.

        Returns:
            Reaction SMARTS for requested template.
        """
        if self.use_db:
            return self.TEMPLATE_DB.find_one({'_id': ObjectId(template_id)})
        else:
            return self.templates[self.id_to_index[template_id]]

    def load_from_database(self):
        """Read the template data from the database."""
        if not self.use_db:
            MyLogger.print_and_log('Error: Cannot load from database when use_db=False',
                transformer_loc, level=3)

        if not self.TEMPLATE_DB:
            self.load_databases()

        # Look for all templates in collection
        to_retrieve = [
            '_id', 'reaction_smarts',
            'necessary_reagent', 'count', 
            'intra_only', 'dimer_only', 'idex',
            'references'
        ]
        for document in self.TEMPLATE_DB.find({}, to_retrieve).sort('index', pymongo.ASCENDING):
            if self.load_all:
                template = self.doc_to_template(document)
                if template is not None:
                    self.templates.append(template)
            else:
                _id = document.get('_id')
                if _id:
                    self.templates.append(_id)
        self.num_templates = len(self.templates)

    def get_outcomes(self, *args, **kwargs):
        """Gets outcome of single transformation.

        Performs a one-step transformation given a SMILES string of a
        target molecule by applying each transformation template
        sequentially.
        """
        raise NotImplementedError

    def load_databases(self, timeout=1000):
        """Loads the databases specified by the global config.

        Args:
            retro (bool): Whether to load the retrosynthetic databases.
            chiral (bool, optional): Whether to properly handle chirality.
                (default: {False})
            timeout (int, optional): Timeout in ms to use before determining the
                database server is not available. (default: {15000})
        """

        db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                                'id'], connect=gc.MONGO['connect'],
                                serverSelectionTimeoutMS=timeout)

        db_name = gc.RETRO_TEMPLATES['database']
        collection = gc.RETRO_TEMPLATES['collection']
        self.TEMPLATE_DB = db_client[db_name][collection]

    def apply_one_template(self, *args, **kwargs):
        """Applies a single template to a given molecule.

        Takes a mol object and applies a single template, returning
        a list of precursors or outcomes, depending on whether retro or
        synthetic templates are used
        """
        raise NotImplementedError
