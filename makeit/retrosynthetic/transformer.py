import os
from rdkit import Chem
import makeit.global_config as gc
from makeit.utilities.io.logger import MyLogger
from makeit.interfaces.template_transformer import TemplateTransformer
from makeit.prioritization.templates.relevance import RelevanceTemplatePrioritizer
from makeit.prioritization.precursors.relevanceheuristic import RelevanceHeuristicPrecursorPrioritizer
from makeit.synthetic.evaluation.fast_filter import FastFilterScorer
from rdchiral.main import rdchiralRun
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from bson.objectid import ObjectId
from pymongo.errors import ServerSelectionTimeoutError

retro_transformer_loc = 'retro_transformer'


class RetroTransformer(TemplateTransformer):
    """Defines an object to perform one-step retrosyntheses for a molecule.

    Attributes:
        mincount (int): Minimum number of precedents for an achiral template for
            inclusion in the template library. Only used when retrotransformers
            need to be initialized.
        mincount_chiral (int): Minimum number of precedents for a chiral
            template for inclusion in the template library. Only used when
            retrotransformers need to be initialized. Chiral templates are
            necessarily more specific, so we generally use a lower threshold
            than achiral templates.
        templates (list): Templates to use for transformation.
        celery (bool): Whether or not Celery is being used.
        TEMPLATE_DB (None or MongoDB): Database to load templates from.
        lookup_only (bool): Whether to only lookup templates in the database
            (instead of loading the entire database).
        precursor_prioritizers (dict): Mapping of precursor prioritizer names to
            objects.
        template_prioritizers (dict):Mapping of template prioritizer names to
            objects.
        precursor_prioritizer (None or str): Specifies which precursor
            prioritizer to use.
        template_prioritizer (None or str): Specifies which template prioritizer
            to use.
        num_templates (int): Number of templates loaded.
        fast_filter (FastFilterScorer or None): Fast filter for evaluation.
        load_all (bool): Whether to load all of the templates into memory.
        use_db (bool): Whether to use the database to look up templates.
        cache_size (int): Maximum cache size to use for template cache. Set to 0
            to not use a cache.
    """


    def __init__(
            self, use_db=True, TEMPLATE_DB=None, load_all=gc.PRELOAD_TEMPLATES,
            template_prioritizer='relevance', 
            precursor_prioritizer='relevanceheuristic',
            fast_filter='default'
        ):
        """Initializes RetroTransformer.

        Args:
            TEMPLATE_DB (None or MongoDB, optional): Database to load
                templates from. (default: {None})
            load_all (bool, optional): Whether to load all of the templates into
                memory. (default: {gc.PRELOAD_TEMPLATES})
            use_db (bool, optional): Whether to use the database to look up
                templates. (default: {True})
        """

        self.templates = []
        self.TEMPLATE_DB = TEMPLATE_DB
        self.template_prioritizer = template_prioritizer
        self.precursor_prioritizer = precursor_prioritizer
        self.fast_filter = fast_filter
        
        super(RetroTransformer, self).__init__(load_all=load_all, use_db=use_db)


    def load(self, template_filename=None):
        if template_filename is None:
            template_filename = os.path.join(
                gc.local_db_dumps,
                gc.RETRO_TRANSFORMS_CHIRAL['file_name']
            )
        if self.template_prioritizer == 'relevance':
            MyLogger.print_and_log('Loading template prioritizer for RetroTransformer', retro_transformer_loc)
            self.template_prioritizer_object = RelevanceTemplatePrioritizer()
            self.template_prioritizer_object.load_model()
            self.template_prioritizer = self.template_prioritizer_object.reorder_templates
        
        if self.precursor_prioritizer == 'relevanceheuristic':
            MyLogger.print_and_log('Loading precursor prioritizer for RetroTransformer', retro_transformer_loc)
            self.precursor_prioritizer_object = RelevanceHeuristicPrecursorPrioritizer()
            self.precursor_prioritizer_object.load_model()
            self.precursor_prioritizer = self.precursor_prioritizer_object.reorder_precursors
        
        if self.fast_filter == 'default':
            MyLogger.print_and_log('Loading fast filter for RetroTransformer', retro_transformer_loc)
            self.fast_filter_object = FastFilterScorer()
            self.fast_filter_object.load(gc.FAST_FILTER_MODEL['trained_model_path'])
            self.fast_filter = lambda x, y: self.fast_filter_object.evaluate(x, y)[0][0]['score']

        MyLogger.print_and_log('Loading retro-synthetic transformer', retro_transformer_loc)
        if self.use_db:
            MyLogger.print_and_log('reading from db', retro_transformer_loc)
            try:
                self.load_from_database()
            except ServerSelectionTimeoutError:
                MyLogger.print_and_log('cannot connect to db, reading from file instead', retro_transformer_loc)
                self.use_db = False
                self.load_from_file(template_filename)
        else:
            MyLogger.print_and_log('reading from file', retro_transformer_loc)
            self.load_from_file(template_filename)

        MyLogger.print_and_log('Retrosynthetic transformer has been loaded - using {} templates.'.format(
            self.num_templates), retro_transformer_loc)

    def get_outcomes(
            self, smiles, template_prioritizer=None, 
            precursor_prioritizer=None, fast_filter=None, 
            fast_filter_threshold=0.75, max_num_templates=100, 
            max_cum_prob=0.995,  **kwargs
        ):
        """Performs a one-step retrosynthesis given a SMILES string.

        Applies each transformation template sequentially to given target
        molecule to perform retrosynthesis.

        Args:
            smiles (str): Target SMILES string to find precursors for.
            template_prioritizer (optional, callable): Use to override
                prioritizer created during initialization. This can be 
                any callable function that accepts
                (smiles, templates, max_num_templates, max_cum_prob) as
                arguments and returns a reordered list of templates
                up until max_num_templates or max_cum_prob as well 
                as a list of the scores for each template.
            precursor_prioritizer (optional, callable): Use to override
                prioritizer created during initialization. This can be
                any callable function that reorders a list of precursor
                dictionary objects.
            fast_filter (optional, callable): Use to override fast filter
                created during initialization. This can be any callable 
                function that accepts (reactants, products) smiles strings 
                as arguments and returns a score on the range [0.0, 1.0].
            fast_filter_threshold (float): Fast filter threshold to filter
                bad predictions. 1.0 means use all templates
            **kwargs: Additional kwargs to pass through to prioritizers or to
                handle deprecated options.

        Returns:
             RetroResult: Special object for a retrosynthetic expansion result,
                defined by ./results.py
        """

        if template_prioritizer is None:
            template_prioritizer = self.template_prioritizer
        
        if precursor_prioritizer is None:
            precursor_prioritizer = self.precursor_prioritizer

        if fast_filter == None:
            fast_filter = self.fast_filter

        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        mol = rdchiralReactants(smiles)

        results = []

        templates, scores = template_prioritizer(smiles, self.templates, max_num_templates, max_cum_prob)

        for template, score in zip(templates, scores):
            if not self.load_all:
                if self.use_db:
                    template = self.TEMPLATE_DB.find_one({'_id': template})
                template = self.doc_to_template(template)
                template['score'] = score
            precursors = self.apply_one_template(mol, template)
            for precursor in precursors:
                joined_smiles = '.'.join(precursor['smiles_list'])
                # skip if no transformation happened
                if joined_smiles == smiles:
                    continue 
                precursor['plasuibility'] = fast_filter(smiles, joined_smiles)
                if precursor['plasuibility'] >= fast_filter_threshold:
                    results.append(precursor)
        results = precursor_prioritizer(results)
        return results

    def apply_one_template(self, mol, template):
        results = []

        try:
            outcomes, mapped_outcomes = rdchiralRun(template['rxn'], mol, return_mapped=True)
        except Exception as e:
            return results

        for j, outcome in enumerate(outcomes):
            smiles_list = []
            smiles_list = outcome.split('.')
            if template['intra_only'] and len(smiles_list) > 1:
                # Disallowed intermolecular reaction
                continue
            if template['dimer_only'] and (len(set(smiles_list)) != 1 or len(smiles_list) != 2):
                # Not a dimer
                continue
            reacting_atoms = mapped_outcomes.get(
                '.'.join(smiles_list), ('.'.join(smiles_list), (-1,))
            )
            results.append({
                'smiles_list': sorted(smiles_list),
                'mapped_smiles': reacting_atoms[0],
                'reacting_atoms': reacting_atoms[1],
                'template_id': str(template['_id']),
                'template_score': template['score'],
                'num_examples': template['count'],
                'necessary_reagent': template['necessary_reagent']
            })
        return results