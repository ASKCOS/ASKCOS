from __future__ import print_function

import makeit.global_config as gc
import os, sys
import makeit.utilities.io.pickle as pickle
from pymongo import MongoClient
from pymongo.errors import ServerSelectionTimeoutError

USE_STEREOCHEMISTRY = True
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from functools import partial  # used for passing args to multiprocessing
from makeit.utilities.io.logger import MyLogger
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

    def __init__(self, celery=False, mincount=gc.RETRO_TRANSFORMS_CHIRAL['mincount'],
        mincount_chiral=gc.RETRO_TRANSFORMS_CHIRAL['mincount_chiral'],
        TEMPLATE_DB=None, lookup_only=False, load_all=gc.PRELOAD_TEMPLATES,
        use_db=True, cache_size=0):
        """Initializes RetroTransformer.

        Args:
            celery (bool, optional): Whether or not Celery is being used.
                (default: {False})
            mincount (int, optional): Minimum number of precedents for an
                achiral template for inclusion in the template library. Only
                used when retrotransformers need to be initialized.
                (default: {25})
            mincount_chiral (int, optional): Minimum number of precedents
                for a chiral template for inclusion in the template library.
                Only used when retrotransformers need to be initialized. Chiral
                templates are necessarily more specific, so we generally use a
                lower threshold than achiral templates. (default: {10})
            TEMPLATE_DB (None or MongoDB, optional): Database to load
                templates from. (default: {None})
            lookup_only (bool, optional): Whether to only lookup templates in
                the database (instead of loading the entire database).
                (default: {False})
            load_all (bool, optional): Whether to load all of the templates into
                memory. (default: {gc.PRELOAD_TEMPLATES})
            use_db (bool, optional): Whether to use the database to look up
                templates. (default: {True})
            cache_size (int, optional): Maximum cache size to use for template
                cache. Set to 0 to not use a cache. (default: {0})
        """

        self.mincount = mincount
        if mincount_chiral == -1:
            self.mincount_chiral = mincount
        else:
            self.mincount_chiral = mincount_chiral
        self.templates = []
        self.celery = celery
        self.TEMPLATE_DB = TEMPLATE_DB
        self.lookup_only = lookup_only
        self.precursor_prioritizers = {}
        self.template_prioritizers = {}
        self.precursor_prioritizer = None
        self.template_prioritizer = None
        self.fast_filter = None
        if self.celery:
            # Pre-load fast filter
            self.load_fast_filter()

        super(RetroTransformer, self).__init__(load_all=load_all, use_db=use_db, cache_size=cache_size)

    def load(self, chiral=True, refs=False, rxns=True, efgs=False, rxn_ex=False):
        """Loads templates to finish initializing the transformer.

        Args:
            TEMPLATE_DB (None or MongoDB, optional): MongoDB to load from.
                (default: {None})
            chiral (bool, optional): Whether to pay close attention to
                chirality. (default: {False})
            refs (bool, optional): Whether to also save references
                (Reaxys instance IDs) when loading templates. (default: {False})
            rxns (bool, optional): Whether to actually load reaction SMARTS
                into RDKit reaction objects. (default: {True})
            efgs (bool, optional): Whether to load statistics about DFG
                popularity. [old] (default: {False})
            rxn_ex (bool, optional): Whether to also save a reaction example
                with each template as it is loaded. (default: {False})
        """

        self.chiral = chiral

        MyLogger.print_and_log('Loading retro-synthetic transformer, including all templates with more than {} hits ({} for chiral reactions)'.format(
            self.mincount, self.mincount_chiral), retro_transformer_loc)

        if chiral:
            from makeit.utilities.io.files import get_retrotransformer_chiral_path
            file_path = get_retrotransformer_chiral_path(
                gc.RETRO_TRANSFORMS_CHIRAL['database'],
                gc.RETRO_TRANSFORMS_CHIRAL['collection'],
                self.mincount,
                self.mincount_chiral,
            )
        else:
            from makeit.utilities.io.files import get_retrotransformer_achiral_path
            file_path = get_retrotransformer_achiral_path(
                gc.RETRO_TRANSFORMS['database'],
                gc.RETRO_TRANSFORMS['collection'],
                self.mincount,
            )

        if self.use_db:
            MyLogger.print_and_log('reading from db', retro_transformer_loc)
            try:
                self.load_from_database(True, chiral=chiral, rxns=True, refs=True, efgs=True, rxn_ex=True)
            except ServerSelectionTimeoutError:
                MyLogger.print_and_log('cannot connect to db, reading from file instead', retro_transformer_loc)
                self.use_db = False
                self.load_from_file(True, file_path, chiral=chiral, rxns=rxns, refs=refs, efgs=efgs, rxn_ex=rxn_ex)
        else:
            MyLogger.print_and_log('reading from file', retro_transformer_loc)
            self.load_from_file(True, file_path, chiral=chiral, rxns=rxns, refs=refs, efgs=efgs, rxn_ex=rxn_ex)

        self.reorder()
        MyLogger.print_and_log('Retrosynthetic transformer has been loaded - using {} templates.'.format(
            self.num_templates), retro_transformer_loc)

    def load_fast_filter(self):
        """Initializes and loads FastFilterScorer.

        NOTE: Keras backend must be Theano for fast filter to work.
        """
        self.fast_filter = FastFilterScorer()
        # self.fast_filter.set_keras_backend('theano')
        self.fast_filter.load(model_path=gc.FAST_FILTER_MODEL['trained_model_path'])

    def get_outcomes(self, smiles, mincount, prioritizers, **kwargs):
        """Performs a one-step retrosynthesis given a SMILES string.

        Applies each transformation template sequentially to given target
        molecule to perform retrosynthesis.

        Args:
            smiles (str): Product SMILES string to find precursors for.
            mincount (int): Minimum template popularity.
            prioritizers (2-tuple of (str, str)): Tuple defining the
                precursor_prioritizer and template_prioritizer to use for
                expansion, each as a string.
            **kwargs: Additional kwargs to pass through to prioritizers or to
                handle deprecated options.

        Returns:
             RetroResult: Special object for a retrosynthetic expansion result,
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
                    filter_flag, filter_score = self.fast_filter.filter_with_threshold(reactant_smiles, smiles, filter_threshold)
                    if filter_flag:
                        precursor.plausibility = filter_score
                        result.add_precursor(precursor, self.precursor_prioritizer, **kwargs)
                else:
                    result.add_precursor(precursor, self.precursor_prioritizer, **kwargs)
        return result

    def apply_one_template_by_idx(self, _id, smiles, template_idx, calculate_next_probs=True, **kwargs):
        """Takes a SMILES and applies the template with given index.

        This is useful in the MCTS code.

        Args:
            _id (int): Not used; passed through to output.
            smiles (str): SMILES of molecule to apply template to.
            template_idx (int): Index of template to be used.
            calculate_next_probs (bool, optional): Whether to calculate template
                relevance probabilities for precursors (default: {True})
            **kwargs: Additional optional arguments.

        Returns:
            list of 5-tuples of (int, str, int, list, float): Result of
                applying given template to the molecule including the template
                relevance probabilities of all resulting precursors when
                calculate_next_probs is True.
        """
        # QUESTION: Why are these not just optional named arguments?
        apply_fast_filter = kwargs.pop('apply_fast_filter', True)
        filter_threshold = kwargs.pop('filter_threshold', 0.75)
        template_count = kwargs.pop('template_count', 100)
        max_cum_prob = kwargs.pop('max_cum_prob', 0.995)
        if (apply_fast_filter and not self.fast_filter):
            self.load_fast_filter()
        self.get_template_prioritizers(gc.relevance)

        # Define mol to operate on
        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)  # to canonicalize
        if self.chiral:
            mol = rdchiralReactants(smiles)

        all_outcomes = []; seen_reactants = {}; seen_reactant_combos = [];
        if self.load_all:
            template = self.templates[template_idx]
        elif not self.use_db:
            template = self.doc_to_template(self.templates[template_idx])
        elif self.cache_size > 0:
            template = self.template_cache[self.templates[template_idx][0]]
        else:
            # HACK: I'm creating a new mongo client every time to avoid the
            #       forking error...
            db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                                    'id'], connect=gc.MONGO['connect'])

            db_name = gc.RETRO_TRANSFORMS_CHIRAL['database']
            collection = gc.RETRO_TRANSFORMS_CHIRAL['collection']
            self.TEMPLATE_DB = db_client[db_name][collection]
            doc = self.TEMPLATE_DB.find_one({'_id': ObjectId(self.templates[template_idx][0])})
            template = self.doc_to_template(doc)
        for smiles_list in self.apply_one_template_smilesonly(mol, smiles, template):
            # Avoid duplicate outcomes (e.g., by symmetry)
            reactant_smiles = '.'.join(smiles_list)
            if reactant_smiles in seen_reactant_combos:
                continue
            seen_reactant_combos.append(reactant_smiles)

            # Should we add this to the results?
            filter_score = 1.0
            if apply_fast_filter:
                filter_flag, filter_score = self.fast_filter.filter_with_threshold(reactant_smiles, smiles, filter_threshold)
                if not filter_flag:
                    continue

            # Should we calculate template relevance scores for each precursor?
            reactants = []
            if calculate_next_probs:
                for reactant_smi in smiles_list:
                    if reactant_smi not in seen_reactants:
                        probs, indeces = self.template_prioritizer.get_topk_from_smi(reactant_smi, k=template_count)
                        # Truncate based on max_cum_prob?
                        truncate_to = np.argwhere(np.cumsum(probs) >= max_cum_prob)
                        if len(truncate_to):
                            truncate_to = truncate_to[0][0] + 1 # Truncate based on max_cum_prob?
                        else:
                            truncate_to = template_count
                        value = 1 # current value assigned to precursor (note: may replace with real value function)
                        # Save to dict
                        seen_reactants[reactant_smi] = (reactant_smi, probs[:truncate_to], indeces[:truncate_to], value)
                    reactants.append(seen_reactants[reactant_smi])

                all_outcomes.append((_id, smiles, template_idx, reactants, filter_score))

            else:
                all_outcomes.append((_id, smiles, template_idx, smiles_list, filter_score))

        if not all_outcomes:
            all_outcomes.append((_id, smiles, template_idx, [], 0.0)) # dummy outcome

        return all_outcomes


    def apply_one_template_smilesonly(self, react_mol, smiles, template, **kwargs):
        """Takes a rdchiralReactants object and applies a single template.

        Only yields SMILES of precursors.

        Args:
            react_mol (rdchiralReactants): Initialized reactant object using
                RDChiral helper package; is the target compound to find
                precursors for.
            smiles (str): Product SMILES (no atom mapping).
            template (dict): Template to be applied, containing an initialized
                rdchiralReaction object as its 'rxn' field.
            smiles_list_only (bool): Whether we only care about the list of
                reactant smiles strings.
            **kwargs: Additional kwargs to accept deprecated options.

        Yields:
            list of str: Results of applying template in the form of a list
                of SMILES strings.
        """
        if template is not None:
            try:
                if self.chiral:
                    outcomes = rdchiralRun(template['rxn'], react_mol)
                else:
                    outcomes = template['rxn'].RunReactants([react_mol])
                results = []
                for j, outcome in enumerate(outcomes):
                    smiles_list = []
                    # Output of rdchiral is (a list of) smiles of the products.
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
                    yield smiles_list
            except Exception as e:
                pass

    def apply_one_template(self, react_mol, smiles, template, **kwargs):
        """Takes a rdchiralReactants object and applies a single template.

        Arguments:
            react_mol (rdchiralReactants): Initialized reactant object using
                RDChiral helper package; is the target compound to find
                precursors for.
            smiles (str): Product SMILES (no atom mapping).
            template (dict): Template to be applied, containing an initialized
                rdchiralReaction object as its 'rxn' field.
            smiles_list_only (bool): Whether we only care about the list of
                reactant smiles strings.
            **kwargs: Additional kwargs to accept deprecated options.

        Returns:
            list: RetroPrecursor objects resulting from applying this one
                template.
        """
        if template is None:
            return []
        try:
            if self.chiral:
                outcomes, mapped_outcomes = rdchiralRun(template['rxn'], react_mol, return_mapped=True)
            else:
                outcomes = template['rxn'].RunReactants([react_mol])
        except Exception as e:
            return []

        results = []
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
                    #cannot have mapped outcomes when not using rdchiral
                    mapped_outcomes = {x:(None,None) for x in smiles_list}
                except Exception as e:
                    print(e) # fail quietly
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

            #Mapped outcomes is {clean_smiles: (mapped_smiles, reacting atoms)}
            reacting_atoms = mapped_outcomes.get('.'.join(smiles_list))

            precursor = RetroPrecursor(
                smiles_list=sorted(smiles_list),
                mapped_smiles=reacting_atoms[0],
                reacting_atoms=reacting_atoms[1],
                template_id=str(template['_id']),
                template_score=template['score'],
                num_examples=template['count'],
                necessary_reagent=template['necessary_reagent']
            )

            results.append(precursor)

        return results

    def top_templates(self, target, **kwargs):
        """Generates only top templates.

        Args:
            target (str): SMILES string of target product.
            **kwargs: Additional options to pass template_prioritizer.

        Yields:
            dict: Single templates in order of decreasing priority.
        """
        for template in self.template_prioritizer.get_priority((self.templates, target),
            TEMPLATE_DB=self.TEMPLATE_DB, mincount=self.mincount,
            mincount_chiral=self.mincount_chiral, chiral=True,
            load_all=self.load_all, use_db=self.use_db,
            template_cache=self.template_cache, **kwargs):

            if not template['chiral'] and template['count'] < self.mincount:
                pass
            elif template['chiral'] and template['count'] < self.mincount_chiral:
                pass
            else:
                yield template

if __name__ == '__main__':

    MyLogger.initialize_logFile()
    t = RetroTransformer()
    t.load(chiral=True, refs=False, rxns=True)

    #Test using a chiral molecule
    outcomes = t.get_outcomes('CCOC(=O)[C@H]1C[C@@H](C(=O)N2[C@@H](c3ccccc3)CC[C@@H]2c2ccccc2)[C@@H](c2ccccc2)N1', \
        100, (gc.relevanceheuristic, gc.relevance))
    precursors = outcomes.precursors

    print([precursor.smiles_list for precursor in precursors])
    print([precursor.reacting_atoms for precursor in precursors])

    #Test using a molecule that give many precursors
    outcomes = t.get_outcomes('CN(C)CCOC(c1ccccc1)c2ccccc2', \
        100, (gc.relevanceheuristic, gc.relevance))
    precursors = outcomes.precursors

    print([precursor.smiles_list for precursor in precursors])
    print([precursor.reacting_atoms for precursor in precursors])
    print([precursor.mapped_smiles for precursor in precursors])

    #test with one template
    outcomes = t.apply_one_template_by_idx(1, 'CCOC(=O)[C@H]1C[C@@H](C(=O)N2[C@@H](c3ccccc3)CC[C@@H]2c2ccccc2)[C@@H](c2ccccc2)N1', 109659)
    print(outcomes)

    # import matplotlib.pyplot as plt
    # import time
    #
    # preload_start = time.time()
    # preload_rt = RetroTransformer(load_all=True)
    # preload_rt.load()
    # preload_elapsed = time.time()-preload_start
    # print('Preload took {} seconds to load.'.format(round(preload_elapsed, 3)))
    # preload_rt.load_fast_filter()
    #
    # db_start = time.time()
    # db_rt = RetroTransformer(load_all=False)
    # db_rt.load()
    # db_elapsed = time.time()-db_start
    # print('DB took {} seconds to load.'.format(round(db_elapsed, 3)))
    # db_rt.load_fast_filter()
    #
    # counts = [10, 100, 1000, 10000, 100000, sys.maxsize]
    # labels = ['10', '100', '1000', '10000', '100000', str(len(db_rt.templates))]
    # smiles = ['CN(C)CCOC(c1ccccc1)c2ccccc2', 'OC(Cn1cncn1)(Cn2cncn2)c3ccc(F)cc3F', 'Cc1ccnc2N(C3CC3)c4ncccc4C(=O)Nc12', 'CN1C2CCC1CC(C2)OC(=O)C(CO)c3ccccc3', 'CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13']
    #
    # counts.reverse()
    # labels.reverse()

    # for p in (1.0, 0.999, 0.997, 0.995):
    #     p_avg = []
    #     p_std = []
    #     db_avg = []
    #     db_std = []
    #     for count in counts:
    #         preload_times = []
    #         db_times = []
    #         for s in smiles:
    #             p_start = time.time()
    #             p_out = preload_rt.get_outcomes(s, 0, (gc.relevanceheuristic, gc.relevance), template_count=count, max_cum_prob=p)
    #             p_elapsed = time.time() - p_start
    #             preload_times.append(p_elapsed)
    #
    #             db_start = time.time()
    #             db_out = db_rt.get_outcomes(s, 0, (gc.relevanceheuristic, gc.relevance), template_count=count, max_cum_prob=p)
    #             db_elapsed = time.time() - db_start
    #             db_times.append(db_elapsed)
    #             print(s)
    #
    #         p_avg.append(np.mean(preload_times))
    #         p_std.append(np.std(preload_times))
    #         db_avg.append(np.mean(db_times))
    #         db_std.append(np.std(db_times))
    #         print(count)
    #
    #     import pickle
    #     with open('{}_p_avg_std_db.pkl'.format(p), 'wb') as f:
    #         pickle.dump((p_avg, p_std, db_avg, db_std), f)
    #     x = np.arange(len(p_avg))
    #     width = 0.35
    #     p1 = plt.bar(x - width/2, p_avg, width, yerr=p_std)
    #     p2 = plt.bar(x + width/2, db_avg, width, yerr=db_std)
    #     plt.ylabel('Elapsed Time (s)')
    #     plt.ylim(bottom=0)
    #     plt.title('Time for Preloading vs DB lookup (p={})'.format(p))
    #     plt.xticks(x, labels)
    #     plt.xlabel('# Templates')
    #     plt.legend((p1[0], p2[0]), ('Preload', 'DB'))
    #     plt.savefig('{}_times.png'.format(p))
    #     plt.close()
    #     print(p)
