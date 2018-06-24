
import makeit.global_config as gc
from multiprocessing import Process, Manager, Queue
import sys 
if sys.version_info[0] < 3:
    import Queue as VanillaQueue
else:
    import queue as VanillaQueue
import time
import sys
from collections import defaultdict
import rdkit.Chem as Chem
from makeit.retrosynthetic.transformer import RetroTransformer
from makeit.utilities.buyable.pricer import Pricer
from makeit.utilities.io.logger import MyLogger
from makeit.utilities.io import model_loader
from makeit.utilities.formats import chem_dict, rxn_dict
import askcos_site.askcos_celery.treebuilder.tb_worker as tb_worker
import askcos_site.askcos_celery.treebuilder.tb_c_worker as tb_c_worker
treebuilder_loc = 'tree_builder'


class TreeBuilder:

    def __init__(self, retroTransformer=None, pricer=None, max_branching=20, max_depth=3, expansion_time=240,
                 celery=False, nproc=1, mincount=25, chiral=True, mincount_chiral=10,
                 template_prioritization=gc.relevance, precursor_prioritization=gc.relevanceheuristic,
                 chemhistorian=None):
        """Class for retrosynthetic tree expansion using a depth-first search

        Initialization of an object of the TreeBuilder class sets default values
        for various settings and loads transformers as needed (i.e., based on 
        whether Celery is being used or not). Most settings are overridden
        by the get_buyable_paths method anyway.

        Keyword Arguments:
            retroTransformer {None or RetroTransformer} -- RetroTransformer object
                to be used for expansion when *not* using Celery. If none, 
                will be initialized using the model_loader.load_Retro_Transformer
                function (default: {None})
            pricer {Pricer} -- Pricer object to be used for checking stop criteria
                (buyability). If none, will be initialized using default settings
                from the global configuration (default: {None})
            max_branching {number} -- Maximum number of precursor suggestions to
                add to the tree at each expansion (default: {20})
            max_depth {number} -- Maximum number of reactions to allow before
                stopping the recursive expansion down one branch (default: {3})
            expansion_time {number} -- Time (in seconds) to allow for expansion
                before searching the generated tree for buyable pathways (default: {240})
            celery {bool} -- Whether or not Celery is being used. If True, then 
                the TreeBuilder relies on reservable retrotransformer workers
                initialized separately. If False, then retrotransformer workers
                will be spun up using multiprocessing (default: {False})
            nproc {number} -- Number of retrotransformer processes to fork for
                faster expansion (default: {1})
            mincount {number} -- Minimum number of precedents for an achiral template
                for inclusion in the template library. Only used when retrotransformers
                need to be initialized (default: {25})
            mincount_chiral {number} -- Minimum number of precedents for a chiral template
                for inclusion in the template library. Only used when retrotransformers
                need to be initialized. Chiral templates are necessarily more specific,
                so we generally use a lower threshold than achiral templates (default: {10})
            chiral {bool} -- Whether or not to pay close attention to chirality. When 
                False, even achiral templates can lead to accidental inversion of
                chirality in non-reacting parts of the molecule. It is highly 
                recommended to keep this as True (default: {True})
            template_prioritization {string} -- Strategy used for template
                prioritization, as a string. There are a limited number of available
                options - consult the global configuration file for info (default: {gc.popularity})
            precursor_prioritization {string} -- Strategy used for precursor
                prioritization, as a string. There are a limited number of available
                options - consult the global configuration file for info (default: {gc.heuristic})
        """
        # General parameters
        self.celery = celery
        self.max_branching = max_branching
        self.mincount = mincount
        self.mincount_chiral = mincount_chiral
        self.max_depth = max_depth
        self.expansion_time = expansion_time
        self.template_prioritization = template_prioritization
        self.precursor_prioritization = precursor_prioritization
        self.nproc = nproc
        self.chiral = chiral
        self.max_cum_template_prob = 1

        if pricer:
            self.pricer = pricer
        else:
            self.pricer = Pricer()
            self.pricer.load()

        self.chemhistorian = chemhistorian
        if chemhistorian is None:
            from makeit.utilities.historian.chemicals import ChemHistorian 
            self.chemhistorian = ChemHistorian()
            self.chemhistorian.load_from_file(refs=False, compressed=True)

        self.reset()

        # When not using Celery, need to ensure retroTransformer initialized
        if not self.celery:
            if retroTransformer:
                self.retroTransformer = retroTransformer
            else:
                self.retroTransformer = model_loader.load_Retro_Transformer(mincount=self.mincount,
                                                                            mincount_chiral=self.mincount_chiral,
                                                                            chiral=self.chiral)

        # Define method to check if all results processed
        if self.celery:
            def waiting_for_results():
                # update
                time.sleep(1)
                return self.pending_results != [] or self.is_ready != []
        else:
            def waiting_for_results():
                waiting = [expansion_queue.empty()
                           for expansion_queue in self.expansion_queues]
                waiting.append(self.results_queue.empty())
                waiting += self.idle

                return (not all(waiting))

        self.waiting_for_results = waiting_for_results

        # Define method to get a processed result.
        if self.celery:
            def get_ready_result():
                # Update which processes are ready
                self.is_ready = [i for (i, res) in enumerate(
                    self.pending_results) if res.ready()]

                for i in self.is_ready:
                    (smiles, precursors) = self.pending_results[
                        i].get(timeout=0.2)
                    self.pending_results[i].forget()
                    _id = self.chem_to_id[smiles]
                    yield (_id, smiles, precursors)
                self.pending_results = [res for (i, res) in enumerate(
                    self.pending_results) if i not in self.is_ready]
        else:
            def get_ready_result():
                while not self.results_queue.empty():
                    yield self.results_queue.get(0.2)

        self.get_ready_result = get_ready_result

        # Define method to start up parallelization
        # note: Celery will reserve an entire preforked pool of workers
        if self.celery:
            def prepare():
                try:
                    if self.chiral:
                        request = tb_c_worker.reserve_worker_pool.delay()
                        self.private_worker_queue = request.get(timeout=10)
                    else:
                        request = tb_worker.reserve_worker_pool.delay()
                        self.private_worker_queue = request.get(timeout=10)
                except Exception as e:
                    request.revoke()
                    raise IOError(
                        'Did not find an available pool of workers! Try again later ({})'.format(e))
        else:
            def prepare():
                MyLogger.print_and_log('Tree builder spinning off {} child processes'.format(
                    self.nproc), treebuilder_loc)
                for i in range(self.nproc):
                    p = Process(target=self.work, args=(i,))
                    self.workers.append(p)
                    p.start()
        self.prepare = prepare

        # Define method to stop working.
        if self.celery:
            def stop():
                if self.pending_results != []:
                    # OPTION 1 - REVOKE TASKS, WHICH GETS SENT TO ALL WORKERS REGARDLESS OF TYPE
                    #[res.revoke() for res in pending_results]
                    # OPTION 2 - DIRECTLY PURGE THE QUEUE (NOTE: HARDCODED FOR
                    # AMQP)
                    import celery.bin.amqp
                    from askcos_site.celery import app
                    amqp = celery.bin.amqp.amqp(app=app)
                    amqp.run('queue.purge', self.private_worker_queue)
                if self.chiral and self.private_worker_queue:
                    released = tb_c_worker.unreserve_worker_pool.apply_async(
                        queue=self.private_worker_queue, retry=True).get()
                elif self.private_worker_queue:
                    released = tb_worker.unreserve_worker_pool.apply_async(
                        queue=self.private_worker_queue, retry=True).get()
                self.running = False
        else:
            def stop():
                if not self.running:
                    return
                self.done.value = 1
                MyLogger.print_and_log(
                    'Terminating tree building process.', treebuilder_loc)

                for p in self.workers:
                    if p and p.is_alive():
                        p.terminate()
                MyLogger.print_and_log(
                    'All tree building processes done.', treebuilder_loc)
                self.running = False
        self.stop = stop

        # Define method to expand a single compound
        if self.celery:
            def expand(smiles, chem_id, depth):
                # Chiral transformation or heuristic prioritization requires
                # same database
                if self.chiral or self.template_prioritization == gc.relevance:
                    self.pending_results.append(tb_c_worker.get_top_precursors.apply_async(
                        args=(smiles, self.template_prioritization,
                              self.precursor_prioritization),
                        kwargs={'mincount': self.mincount,
                                'max_branching': self.max_branching,
                                'template_count': self.template_count,
                                'mode': self.precursor_score_mode,
                                'max_cum_prob': self.max_cum_template_prob,
                                'apply_fast_filter': self.apply_fast_filter,
                                'filter_threshold': self.filter_threshold},
                        # Prioritize higher depths: Depth first search.
                        priority=int(depth),
                        queue=self.private_worker_queue,
                    ))
                else:
                    self.pending_results.append(tb_worker.get_top_precursors.apply_async(
                        args=(smiles, self.template_prioritization,
                              self.precursor_prioritization),
                        kwargs={'mincount': self.mincount,
                                'max_branching': self.max_branching,
                                'template_count': self.template_count,
                                'mode': self.precursor_score_mode,
                                'max_cum_prob': self.max_cum_template_prob,
                                'apply_fast_filter': self.apply_fast_filter,
                                'filter_threshold': self.filter_threshold},
                        # Prioritize higher depths: Depth first search.
                        priority=int(depth),
                        queue=self.private_worker_queue,
                    ))
        else:
            def expand(smiles, chem_id, depth):
                #print('Coordinator put {} (ID {}) in queue queue {}'.format(smiles, chem_id, depth))
                self.expansion_queues[depth].put((chem_id, smiles))
        self.expand = expand

        # Define how first target is set.
        if self.celery:
            def set_initial_target(smiles):
                self.expand(smiles, 1, 0)
        else:
            def set_initial_target(smiles):
                self.expansion_queues[-1].put((1, smiles))
                print(self.expansion_queues)
                print('Put something on expansion queue')
                while self.results_queue.empty():
                    time.sleep(0.25)
                    #print('Waiting for first result in treebuilder...')
        self.set_initial_target = set_initial_target

    def reset(self):
        """Clears saved state and resets counters

        Called after initialization and after getting buyable pathways
        to free up memory and - in the case of Celery - be prepared to
        handle another request without having results carry over between
        different tree building requests.
        """
        if self.celery:
            # general parameters in celery format
            self.tree_dict = {}
            self.chem_to_id = {}
            self.buyable_leaves = set()
            self.current_id = 2
            self.is_ready = []
            # specifically for celery
            self.pending_results = []
            self.private_worker_queue = None
        else:
            # general parameters in python multiprocessing format
            self.manager = Manager()
            self.tree_dict = self.manager.dict()
            self.chem_to_id = self.manager.dict()
            self.buyable_leaves = self.manager.list()
            self.current_id = self.manager.Value('i', 2)

            # specificly for python multiprocessing
            self.done = self.manager.Value('i', 0)
            self.paused = self.manager.Value('i', 0)
            # Keep track of idle workers
            self.idle = self.manager.list()
            self.results_queue = Queue()
            self.workers = []
            self.coordinator = None
            self.running = False

    def get_children(self, precursors):
        """Reformats result of an expansion for tree buidler

        Arguments:
            precursors {list of dicts} -- precursor information as a list of
                dictionaries, where each dictionary contains information about
                the precursor identity and auxiliary information like
                whether the template requires heavy atoms to be contributed
                by reagents.

        Returns:
            [list of (dict, list)] -- precursor information reformatted into
                a list of 2-tuples, where the 2-tuple consists of a dictionary
                containing basically the same information as the precursor 
                dictionaries, as well as the smiles_split of the precursor,
                which is a list of reactant SMILES strings
        """
        children = []
        for precursor in precursors:
            children.append((
                {
                    'tforms': precursor['tforms'],
                    'template': precursor['tforms'][0],
                    'template_score': precursor['template_score'],
                    'necessary_reagent': precursor['necessary_reagent'],
                    'num_examples': precursor['num_examples'],
                    'score': precursor['score'],
                    'plausibility': precursor['plausibility']
                },
                precursor['smiles_split']
            ))

        return children

    def add_children(self, children, smiles, unique_id):
        """Add results of one expansion to the tree dictionary

        Arguments:
            children {list of (dict, list)} -- list of candidate disconnections
                for the target SMILES, formatted as returned by get_children()
            smiles {string} -- SMILES string of product (target) molecule
            unique_id {integer >= 1} -- ID of product (target) molecule in the
                tree dictionary
        """

        parent_chem_doc = self.tree_dict[unique_id]  # copy to overwrite later
        parent_chem_prod_of = parent_chem_doc['prod_of']
        # Assign unique number
        for (rxn, mols) in children:

            # Add option to leave out blacklisted reactions.
            rxn_smiles = '.'.join(sorted(mols)) + '>>' + smiles
            if rxn_smiles in self.known_bad_reactions:
                continue

            # What should be excluded?
            skip_this = False
            for mol in mols:
                # Exclude banned molecules too
                if mol in self.forbidden_molecules:
                    skip_this = True
                # Exclude reactions where the reactant is the target
                if mol == self.tree_dict[1]['smiles']:
                    skip_this = True
            if skip_this:
                continue

            # depending on whether current_id was given as 'Manager.Value' type
            # or 'Integer':
            if self.celery:
                rxn_id = self.current_id
                self.current_id += 1
            else:
                rxn_id = self.current_id.value
                # this is only okay because there is/should be only ONE
                # treebuilder
                self.current_id.value += 1
            # For the parent molecule, record child reactions
            parent_chem_prod_of.append(rxn_id)

            # For the reaction, keep track of children IDs
            chem_ids = []
            for mol in mols:

                # New chemical?
                if mol not in self.chem_to_id:

                    try:
                        chem_id = self.current_id.value
                        # this is only okay because there is/should be only ONE
                        # treebuilder
                        self.current_id.value += 1
                    except AttributeError:
                        chem_id = self.current_id
                        self.current_id += 1

                    # Check if buyable
                    ppg = self.pricer.lookup_smiles(mol, alreadyCanonical=True)
                    hist = self.chemhistorian.lookup_smiles(mol, alreadyCanonical=True)

                    self.tree_dict[chem_id] = {
                        'smiles': mol,
                        'prod_of': [],
                        'rct_of': [rxn_id],
                        'depth': parent_chem_doc['depth'] + 1,
                        'ppg': ppg,
                        'as_reactant': hist['as_reactant'],
                        'as_product': hist['as_product'],
                    }
                    self.chem_to_id[mol] = chem_id

                    # Check stop criterion
                    if self.is_a_leaf_node(mol, ppg, hist):
                        # print('{} is a leaf!'.format(mol))
                        if self.celery:
                            self.buyable_leaves.add(chem_id)
                        else:
                            self.buyable_leaves.append(chem_id)

                    else:
                        # Add to queue to get expanded
                        if parent_chem_doc['depth'] >= self.max_depth - 1:
                            if gc.DEBUG:
                                MyLogger.print_and_log('Reached maximum depth, so will not expand around {}'.format(
                                    tree_dict[chem_id]), treebuilder_loc)
                        else:
                            self.expand(mol, chem_id, parent_chem_doc['depth'])

                else:
                    chem_id = self.chem_to_id[mol]

                    # Overwrite this chemical node to record it is a reactant
                    # of this rxn
                    chem_doc = self.tree_dict[chem_id]
                    chem_doc['rct_of'] += [rxn_id]
                    self.tree_dict[chem_id] = chem_doc

                # Save ID
                chem_ids.append(chem_id)

            # Record by overwriting the whole dict value
            rxn['rcts'] = chem_ids
            rxn['prod'] = unique_id
            rxn['depth'] = parent_chem_doc['depth'] + 0.5
            self.tree_dict[rxn_id] = rxn

        # Overwrite dictionary entry for the parent
        parent_chem_doc['prod_of'] = parent_chem_prod_of
        self.tree_dict[unique_id] = parent_chem_doc

    def work(self, i):
        """Work function for retroTransformer processes

        Only used when Celery is false and multiprocessing is used instead. Will
        constantly look for molecules on expansion queues to expand (in a DFS) 
        and add results to the results queue.

        Arguments:
            i {integer >= 0} -- index assigned to the worker, used to assign
                idle status to the shared list self.idle[i]
        """
        while True:
            # If done, stop
            if self.done.value:
                MyLogger.print_and_log(
                    'Worker {} saw done signal, terminating'.format(i), treebuilder_loc)
                break
            # If paused, wait and check again
            if self.paused.value:
                #print('Worker {} saw pause signal, sleeping for 1 second'.format(i))
                time.sleep(1)
                continue
            # Grab something off the queue
            for j in range(len(self.expansion_queues))[::-1]:
                try:
                    (_id, smiles) = self.expansion_queues[j].get(timeout=0.1)  # short timeout
                    self.idle[i] = False
                    # print('Worker {} grabbed {} (ID {}) to expand from queue {}'.format(i, smiles, _id, j))
                    result = self.retroTransformer.get_outcomes(smiles, self.mincount, (self.precursor_prioritization,
                                                                                        self.template_prioritization),
                                                                template_count=self.template_count,
                                                                mode=self.precursor_score_mode,
                                                                max_cum_prob=self.max_cum_template_prob,
                                                                apply_fast_filter=self.apply_fast_filter,
                                                                filter_threshold=self.filter_threshold
                                                                )

                    precursors = result.return_top(n=self.max_branching)

                    self.results_queue.put((_id, smiles, precursors))
                    

                except VanillaQueue.Empty:
                    #print('Queue {} empty for worker {}'.format(j, i))
                    pass
                except Exception as e:
                    sys.stdout.write(e)
                    sys.stdout.flush()
            time.sleep(0.01)
            self.idle[i] = True

    def coordinate(self):
        """Run the expansion up to a specified self.expansion_time (seconds)
        """
        start_time = time.time()
        elapsed_time = time.time() - start_time
        next = 1
        while (elapsed_time < self.expansion_time) and self.waiting_for_results():
            if (int(elapsed_time)/10 == next):
                next += 1
                MyLogger.print_and_log(
                    'Worked for {}/{} s'.format(int(elapsed_time*10)/10.0, self.expansion_time), treebuilder_loc)
            try:
                for (_id, smiles, precursors) in self.get_ready_result():
                    children = self.get_children(precursors)
                    self.add_children(children, smiles, _id)
                elapsed_time = time.time() - start_time
            except Exception as e:
                elapsed_time = time.time() - start_time
                print('##ERROR#: {}'.format(e))

    def build_tree(self, target):
        """Recursively build out the synthesis tree

        Arguments:
            target {string} -- SMILES of target molecule
        """
        self.running = True
        if self.celery:
            from celery.result import allow_join_result
        else:
            from makeit.utilities.with_dummy import with_dummy as allow_join_result
        with allow_join_result():
            try:
                hist = self.chemhistorian.lookup_smiles(target)
                self.tree_dict[1] = {
                    'smiles': target,
                    'prod_of': [],
                    'rct_of': [],
                    'depth': 0,
                    'ppg': self.pricer.lookup_smiles(target),
                    'as_reactant': hist['as_reactant'],
                    'as_product': hist['as_product'],
                }

                if self.is_a_leaf_node(target, self.tree_dict[1]['ppg'], hist):
                    if self.celery:
                        self.buyable_leaves.add(1)
                    else:
                        self.buyable_leaves.append(1)

                self.chem_to_id[target] = 1

                self.prepare()
                self.set_initial_target(target)
                self.coordinate()

            finally:  # make sure stop is graceful
                self.stop()

    def tree_status(self):
        """Summarize size of tree after expansion

        Returns:
            num_chemicals {int} -- number of chemical nodes in the tree
            num_reactions {int} -- number of reaction nodes in the tree
            at_depth {dict} -- dictionary containing counts at each integer
                depth (chemicals) and half-integer depth (reactions)
        """

        num_chemicals = 0
        num_reactions = 0
        at_depth = {}
        for _id in self.tree_dict.keys():
            depth = self.tree_dict[_id]['depth']
            if depth % 1 == 0:
                num_chemicals += 1
            else:
                num_reactions += 1
            if depth not in at_depth:
                at_depth[depth] = 1
            else:
                at_depth[depth] += 1
        return (num_chemicals, num_reactions, at_depth)

    def get_buyable_paths(self, target, max_depth=3, max_branching=25, expansion_time=240, template_prioritization=gc.relevance,
                          precursor_prioritization=gc.heuristic, nproc=1, mincount=25, chiral=True, mincount_chiral=10, max_trees=25, max_ppg=1e10,
                          known_bad_reactions=[], forbidden_molecules=[], template_count=100, precursor_score_mode=gc.max,
                          max_cum_template_prob=1, max_natom_dict=defaultdict(lambda: 1e9, {'logic': None}),
                          min_chemical_history_dict={'as_reactant':1e9, 'as_product':1e9,'logic':None},
                          apply_fast_filter= False, filter_threshold=0.5):
        """Get viable synthesis trees using an iterative deepening depth-first search

        [description]

        Arguments:
            target {[type]} -- [description]

        Keyword Arguments:
            max_depth {number} -- Maximum number of reactions to allow before
                stopping the recursive expansion down one branch (default: {3})
            max_branching {number} -- Maximum number of precursor suggestions to
                add to the tree at each expansion (default: {25})
            expansion_time {number} -- Time (in seconds) to allow for expansion
                before searching the generated tree for buyable pathways (default: {240})
            nproc {number} -- Number of retrotransformer processes to fork for
                faster expansion (default: {1})
            mincount {number} -- Minimum number of precedents for an achiral template
                for inclusion in the template library. Only used when retrotransformers
                need to be initialized (default: {25})
            mincount_chiral {number} -- Minimum number of precedents for a chiral template
                for inclusion in the template library. Only used when retrotransformers
                need to be initialized. Chiral templates are necessarily more specific,
                so we generally use a lower threshold than achiral templates (default: {10})
            chiral {bool} -- Whether or not to pay close attention to chirality. When 
                False, even achiral templates can lead to accidental inversion of
                chirality in non-reacting parts of the molecule. It is highly 
                recommended to keep this as True (default: {True})
            template_prioritization {string} -- Strategy used for template
                prioritization, as a string. There are a limited number of available
                options - consult the global configuration file for info (default: {gc.relevance})
            precursor_prioritization {string} -- Strategy used for precursor
                prioritization, as a string. There are a limited number of available
                options - consult the global configuration file for info (default: {gc.heuristic})
            max_trees {number} -- Maximum number of buyable trees to return. Does
                not affect expansion time (default: {25})
            max_ppg {number} -- Maximum price ($/g) for a chemical to be considered
                buyable, and thus potentially usable as a leaf node (default: {1e10})
            known_bad_reactions {list} -- Reactions to forbid during expansion, 
                represented as list of reaction SMILES strings. Each reaction SMILES
                must be canonicalized, have atom mapping removed, and have its
                reactant fragments be sorted. Forbidden reactions are checked 
                when processing children returned by the RetroTransformer (default: {[]})
            forbidden_molecules {list} -- Molecules to forbid during expansion, 
                represented as a list of SMILES strings. Each string must be
                canonicalized without atom mapping. Forbidden molecules will not
                be allowed as intermediates or leaf nodes (default: {[]})
            template_count {number} -- Maximum number of templates to apply at
                each expansion (default: {100})
            precursor_score_mode {string} -- Mode to use for precursor scoring
                when using the SCScore prioritizer and multiple reactant
                fragments must be scored together (default: {gc.max})
            max_cum_template_prob {number} -- Maximum cumulative template
                probability (i.e., relevance score), used as part of the
                relevance template_prioritizer (default: {1})
            max_natom_dict {dict} -- Dictionary defining a potential chemical
                property stopping criterion based on the number of atoms of
                C, N, O, and H in a molecule. The 'logic' keyword of the dict
                refers to how that maximum number of atom information is combined
                with the requirement that chemicals be cheaper than max_ppg
                (default: {defaultdict(lambda: 1e9, {'logic': None})})
            min_chemical_history_dict {dict} -- Dictionary defining a potential
                chemical stopping criterion based on the number of times a
                molecule has been seen previously. Always uses logical OR

        Returns:
            tree_status -- result of tree_status()
            trees -- list of dictionaries, where each dictionary defines a
                synthetic route
        """

        self.mincount = mincount
        self.mincount_chiral = mincount_chiral
        self.max_depth = max_depth
        self.max_branching = max_branching
        self.expansion_time = expansion_time
        self.template_prioritization = template_prioritization
        self.precursor_prioritization = precursor_prioritization
        self.precursor_score_mode = precursor_score_mode
        self.nproc = nproc
        self.template_count = template_count
        self.max_cum_template_prob = max_cum_template_prob
        self.max_ppg = max_ppg
        self.apply_fast_filter= apply_fast_filter
        self.filter_threshold = filter_threshold

        MyLogger.print_and_log('Starting to expand {} using max_natom_dict {}, min_history {}'.format(
            target, max_natom_dict, min_chemical_history_dict), treebuilder_loc, level=1)

        if min_chemical_history_dict['logic'] not in [None, 'none'] and \
                self.chemhistorian is None:
            from makeit.utilities.historian.chemicals import ChemHistorian 
            self.chemhistorian = ChemHistorian()
            self.chemhistorian.load_from_file(refs=False, compressed=True)
            MyLogger.print_and_log('Loaded compressed chemhistorian from file', treebuilder_loc, level=1)

        # Define stop criterion
        def is_buyable(ppg):
            return ppg and (ppg <= self.max_ppg)
        def is_small_enough(smiles):
            # Get structural properties
            natom_dict = defaultdict(lambda: 0)
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return False
            for a in mol.GetAtoms():
                natom_dict[a.GetSymbol()] += 1
            natom_dict['H'] = sum(a.GetTotalNumHs() for a in mol.GetAtoms())
            max_natom_satisfied = all(natom_dict[k] <= v for (
                k, v) in max_natom_dict.items() if k != 'logic')
            return max_natom_satisfied
        def is_popular_enough(hist):
            return hist['as_reactant'] >= min_chemical_history_dict['as_reactant'] or \
                    hist['as_product'] >= min_chemical_history_dict['as_product']
        
        if min_chemical_history_dict['logic'] in [None, 'none']:
            if max_natom_dict['logic'] in [None, 'none']:
                def is_a_leaf_node(smiles, ppg, hist):
                    return is_buyable(ppg)
            elif max_natom_dict['logic'] == 'or':
                def is_a_leaf_node(smiles, ppg, hist):
                    return is_buyable(ppg) or is_small_enough(smiles)
            else:
                def is_a_leaf_node(smiles, ppg, hist):
                    return is_buyable(ppg) and is_small_enough(smiles)
        else:
            if max_natom_dict['logic'] in [None, 'none']:
                def is_a_leaf_node(smiles, ppg, hist):
                    return is_buyable(ppg) or is_popular_enough(hist)
            elif max_natom_dict['logic'] == 'or':
                def is_a_leaf_node(smiles, ppg, hist):
                    return is_buyable(ppg) or is_popular_enough(hist) or is_small_enough(smiles)
            else:
                def is_a_leaf_node(smiles, ppg, hist):
                    return is_popular_enough(hist) or (is_buyable(ppg) and is_small_enough(smiles))
            
        self.is_a_leaf_node = is_a_leaf_node

        # Override: if relevance method is used, chiral database must be used!
        if chiral or template_prioritization == gc.relevance:
            self.chiral = True
        if template_prioritization == gc.relevance and not self.celery:
            if not (self.retroTransformer.mincount == 25
                    and self.retroTransformer.mincount_chiral == 10
                    and self.retroTransformer.chiral):
                MyLogger.print_and_log('When using relevance based template prioritization, chiral template database ' +
                                       'must be used with mincount = 25 and mincount_chiral = 10. Exiting...', treebuilder_loc, level=3)

        self.known_bad_reactions = known_bad_reactions
        self.forbidden_molecules = forbidden_molecules
        self.reset()

        # Initialize multiprocessing queues if necessary
        if not self.celery:
            for i in range(nproc):
                self.idle.append(True)
            if self.max_depth != 1:
                self.expansion_queues = [Queue()
                                         for i in range(self.max_depth - 1)]
            else:
                self.expansion_queues = [Queue()]

        # Generate trees
        self.build_tree(target)

        def IDDFS():
            """Perform an iterative deepening depth-first search to find buyable
            pathways.
                        
            Yields:
                nested dictionaries defining synthesis trees
            """
            for depth in range(self.max_depth+1):
                for path in DLS_chem(1, depth, headNode=True):
                    yield chem_dict(1, children=path, **self.tree_dict[1])

        def DLS_chem(chem_id, depth, headNode=False):
            """Expand at a fixed depth for the current node chem_id.
            
            Arguments:
                chem_id {int >= 1} -- unique ID of the current chemical
                depth {int >= 0} -- current depth of the search
            
            Keyword Arguments:
                headNode {bool} -- Whether this is the first (head) node, in which
                    case it must be expanded even if it is buyable itself (default: {False})

            Yields:
                list -- paths connecting to buyable molecules that are children
                    of the current chemical
            """
            # Copy list so each new branch has separate list.

            if depth <= 0:
                # Not allowing deeper - is this buyable?
                if chem_id in self.buyable_leaves:
                    yield []  # viable node, calling function doesn't need children
            else:
                # Do we need to go deeper?
                if chem_id in self.buyable_leaves and not headNode:
                    yield []  # Nope, this is a viable node
                else:
                    # Try going deeper via DLS_rxn function
                    for rxn_id in self.tree_dict[chem_id]['prod_of']:
                        rxn_smiles = '.'.join(sorted([self.tree_dict[x]['smiles'] for x in self.tree_dict[
                                            rxn_id]['rcts']])) + '>>' + self.tree_dict[chem_id]['smiles']
                        for path in DLS_rxn(rxn_id, depth):
                            yield [rxn_dict(rxn_id, rxn_smiles, children=path, **self.tree_dict[rxn_id])]

        def DLS_rxn(rxn_id, depth):
            """Return children paths starting from a specific rxn_id
            
            Arguments:
                rxn_id {int >= 2} -- unique ID of this reaction in the tree_dict
                depth {int >= 0} -- current depth of the search
            
            Yields:
                list -- paths connecting to buyable molecules that are children
                    of the current reaction
            """

            # Only one reactant? easy!
            if len(self.tree_dict[rxn_id]['rcts']) == 1:
                chem_id = self.tree_dict[rxn_id]['rcts'][0]
                for path in DLS_chem(chem_id, depth-1):
                    yield [
                        chem_dict(chem_id, children=path, **self.tree_dict[chem_id])
                    ]

            # Two reactants? want to capture all combinations of each node's
            # options
            elif len(self.tree_dict[rxn_id]['rcts']) == 2:
                chem_id0 = self.tree_dict[rxn_id]['rcts'][0]
                chem_id1 = self.tree_dict[rxn_id]['rcts'][1]
                for path0 in DLS_chem(chem_id0, depth-1):
                    for path1 in DLS_chem(chem_id1, depth-1):
                        yield [
                            chem_dict(chem_id0, children=path0, **self.tree_dict[chem_id0]),
                            chem_dict(chem_id1, children=path1, **self.tree_dict[chem_id1]),
                        ]

            # Three reactants? This is not elegant...
            elif len(self.tree_dict[rxn_id]['rcts']) == 3:
                chem_id0 = self.tree_dict[rxn_id]['rcts'][0]
                chem_id1 = self.tree_dict[rxn_id]['rcts'][1]
                chem_id2 = self.tree_dict[rxn_id]['rcts'][2]
                for path0 in DLS_chem(chem_id0, depth-1):
                    for path1 in DLS_chem(chem_id1, depth-1):
                        for path2 in DLS_chem(chem_id2, depth-1):
                            yield [
                                chem_dict(chem_id0, children=path0, **self.tree_dict[chem_id0]),
                                chem_dict(chem_id1, children=path1, **self.tree_dict[chem_id1]),
                                chem_dict(chem_id2, children=path2, **self.tree_dict[chem_id2]),
                            ]

            # I am ashamed
            elif len(self.tree_dict[rxn_id]['rcts']) == 4:
                chem_id0 = self.tree_dict[rxn_id]['rcts'][0]
                chem_id1 = self.tree_dict[rxn_id]['rcts'][1]
                chem_id2 = self.tree_dict[rxn_id]['rcts'][2]
                chem_id3 = self.tree_dict[rxn_id]['rcts'][3]
                for path0 in DLS_chem(chem_id0, depth-1):
                    for path1 in DLS_chem(chem_id1, depth-1):
                        for path2 in DLS_chem(chem_id2, depth-1):
                            for path3 in DLS_chem(chem_id3, depth-1):
                                yield [
                                    chem_dict(chem_id0, children=path0, **self.tree_dict[chem_id0]),
                                    chem_dict(chem_id1, children=path1, **self.tree_dict[chem_id1]),
                                    chem_dict(chem_id2, children=path2, **self.tree_dict[chem_id2]),
                                    chem_dict(chem_id3, children=path3, **self.tree_dict[chem_id3]),
                                ]

            else:
                print('Too many reactants! Only have cases 1-4 programmed')
                raise ValueError(
                    'Too many reactants! Only have cases 1-4 programmed')

        # Generate paths and ensure unique
        import hashlib
        import json
        done_trees = set()
        trees = []
        counter = 0
        for tree in IDDFS():
            hashkey = hashlib.sha1(json.dumps(
                tree, sort_keys=True).encode('utf-8')).hexdigest()

            if hashkey in done_trees:
                #print('Found duplicate tree...')
                continue

            done_trees.add(hashkey)
            trees.append(tree)
            counter += 1

            if counter == max_trees:
                MyLogger.print_and_log('Generated {} trees (max_trees met), stopped looking for more...'.format(
                    max_trees), treebuilder_loc)
                break

        tree_status = self.tree_status()
        if self.celery:
            self.reset()  # free up memory, don't hold tree
        return (tree_status, trees)

if __name__ == '__main__':
    MyLogger.initialize_logFile()
    celery = False
    treeBuilder = TreeBuilder(celery=celery, mincount=25, mincount_chiral=10)
    status, paths = treeBuilder.get_buyable_paths('CN1C2CCC1CC(OC(=O)C(CO)c1ccccc1)C2', max_depth=4, template_prioritization=gc.relevance,
                                        precursor_prioritization=gc.relevanceheuristic, nproc=2, expansion_time=60, max_trees=500, max_ppg=10,
                                        max_branching=25, apply_fast_filter=True, filter_threshold=0.75,
                                        min_chemical_history_dict={'as_reactant':5, 'as_product':1, 'logic':'none'})
    print('done')
    print(status)
    print(paths[0])
    print(paths[1])
    print(paths[2])
