import makeit.global_config as gc
from makeit.utilities.io.logger import MyLogger
from rdchiral.main import rdchiralRun
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from bson.objectid import ObjectId

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


    def __init__(self, TEMPLATE_DB=None, load_all=gc.PRELOAD_TEMPLATES):
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
        super(RetroTransformer, self).__init__(load_all=load_all, use_db=use_db)


    def load(self, template_filename=gc.RETRO_TRANSFORMS_CHIRAL['file_name']):
        MyLogger.print_and_log('Loading retro-synthetic transformer', retro_transformer_loc)
        if self.use_db:
            MyLogger.print_and_log('reading from db', retro_transformer_loc)
            try:
                self.load_from_database()
            except ServerSelectionTimeoutError:
                MyLogger.print_and_log('cannot connect to db, reading from file instead', retro_transformer_loc)
                self.use_db = False
                self.load_from_file(True, template_filename)
        else:
            MyLogger.print_and_log('reading from file', retro_transformer_loc)
            self.load_from_file(True, file_path)

        MyLogger.print_and_log('Retrosynthetic transformer has been loaded - using {} templates.'.format(
            self.num_templates), retro_transformer_loc)

    def get_outcomes(self, smiles, template_prioritizer='relevance', precursor_prioritizer='relevanceheuristic', fast_filter='default', fast_filter_threshold=0.75,  **kwargs):
        """Performs a one-step retrosynthesis given a SMILES string.

        Applies each transformation template sequentially to given target
        molecule to perform retrosynthesis.

        Args:
            smiles (str): Target SMILES string to find precursors for.
            template_prioritizer (str or function): method for template 
                prioritization. This can either be 'relevance' or any 
                callable function that accepts (smiles, templates) as
                arguments and returns a reordered list of templates.
            precursor_prioritizer (str or function): method for precursor
                prioritization. This can either be 'relevanceheuristic' or
                any callable function that reorders a list of precursor
                smiles strings.
            fast_filter (str of function): method for filtering precursors.
                This can either be 'default' or any function that accepts
                (target, precursor) smiles strings as arguments and returns
                a score on the range [0.0, 1.0].
            fast_filter_threshold (float): Fast filter threshold to filter
                bad predictions. 1.0 means use all templates
            **kwargs: Additional kwargs to pass through to prioritizers or to
                handle deprecated options.

        Returns:
             RetroResult: Special object for a retrosynthetic expansion result,
                defined by ./results.py
        """

        if template_prioritizer == 'relevance':
            template_prioritizer = #TODO
        
        if precursor_prioritizer == 'relevanceheuristic':
            precursor_prioritizer = #TODO

        if fast_filter == 'default':
            fast_filter = #TODO

        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        mol = rdchiralReactants(smiles)

        results = []

        for template in template_prioritizer(smiles, self.templates):
            for precursor in self.apply_one_template(mol, smiles, template):
                precursor['plasuibility'] = fast_filter(smiles, '.'.join(precursor['smiles_list']))
                if precursor['plasuibility'] >= fast_filter_threshold:
                    results.append(precursor)
        results = precursor_prioritizer(results)
        return results

    def apply_one_template(self, mol, template):
        results = []
        try:
            outcomes, mapped_outcomes = rdchiralRun(template['rxn'], mol, return_mapped=True)
        except Exception as e:
                pass
        for j, outcome in enumerate(outcomes):
            smiles_list = []
            smiles_list = outcome.split('.')
            if template['intra_only'] and len(smiles_list) > 1:
                # Disallowed intermolecular reaction
                continue
            if template['dimer_only'] and (len(set(smiles_list)) != 1 or len(smiles_list) != 2):
                # Not a dimer
                continue
            if '.'.join(smiles_list) == smiles:
                # no transformation
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