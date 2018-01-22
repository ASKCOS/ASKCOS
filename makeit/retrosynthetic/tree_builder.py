
import makeit.global_config as gc
from multiprocessing import Process, Manager, Queue
import Queue as VanillaQueue
import time
import sys
from makeit.retrosynthetic.transformer import RetroTransformer
from makeit.utilities.buyable.pricer import Pricer
from makeit.utilities.io.logging import MyLogger
from makeit.utilities.io import model_loader
from makeit.utilities.formats import chem_dict, rxn_dict
from celery.result import allow_join_result
import askcos_site.askcos_celery.treebuilder.tb_worker as tb_worker
import askcos_site.askcos_celery.treebuilder.tb_c_worker as tb_c_worker
treebuilder_loc = 'tree_builder'


class TreeBuilder:
    '''
    Class for retro-synthetic tree expansion. The expansion is performed based on the retroTransformer that is provided.
    '''

    def __init__(self, retroTransformer=None, pricer=None, max_branching=20, max_depth=3, expansion_time=240,
                 celery=False, nproc=1, mincount=25, chiral=True, template_prioritization=gc.popularity,
                 precursor_prioritization=gc.heuristic, mincount_chiral=10):
        # General parameters
        self.celery = celery
        self.max_branching = max_branching
        self.mincount = mincount
        self.mincount_chiral = mincount_chiral
        self.max_depth = max_depth
        self.expansion_time = expansion_time
        self.template_prioritization = template_prioritization
        MyLogger.print_and_log('Using {} method for template prioritization.'.format(
            self.template_prioritization), treebuilder_loc)
        self.precursor_prioritization = precursor_prioritization
        MyLogger.print_and_log('Using {} method for precursor prioritization.'.format(
            self.precursor_prioritization), treebuilder_loc)
        self.nproc = nproc
        self.chiral = chiral

        if pricer:
            self.pricer = pricer
        else:
            self.pricer = Pricer()

        self.reset()

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
                #Update which processes are ready
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

        # Define method to start up parallelization.
        if self.celery:
            def prepare():
                if self.chiral:
                    self.private_worker_queue = tb_c_worker.reserve_worker_pool.delay().get(timeout=5)
                else:
                    self.private_worker_queue = tb_worker.reserve_worker_pool.delay().get(timeout=5)
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
                if self.chiral:
                    released = tb_c_worker.unreserve_worker_pool.apply_async(
                        queue=self.private_worker_queue, retry=True).get()
                else:
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
                                'mode': self.precursor_score_mode},
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
                                'mode': self.precursor_score_mode},
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
                while self.results_queue.empty():
                    time.sleep(0.25)
        self.set_initial_target = set_initial_target

    def reset(self):
        if self.celery:
            # general parameters in celery format
            self.tree_dict = {}
            self.chem_to_id = {}
            self.buyable_leaves = set()
            self.current_id = 2
            self.is_ready = []
            # specifically for celery
            self.pending_results = []
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
        '''
        Format children for a given set of precursors.
        '''
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
                },
                precursor['smiles_split']
            ))

        return children

    def add_children(self, children, smiles, unique_id):
        '''
        Add the precursors to the dictionary the builder was initiated with. Should be used in the "Coordinator" in
        in parallel processing.
        To ensure correct working of the indexation, only one treebuilder object should be instantiated.
        '''
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

                    self.tree_dict[chem_id] = {
                        'smiles': mol,
                        'prod_of': [],
                        'rct_of': [rxn_id],
                        'depth': parent_chem_doc['depth'] + 1,
                        'ppg': ppg
                    }
                    self.chem_to_id[mol] = chem_id

                    if ppg:
                        #print('{} buyable!'.format(mol))
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
                    (_id, smiles) = self.expansion_queues[
                        j].get(timeout=0.5)  # short timeout
                    self.idle[i] = False
                    #print('Worker {} grabbed {} (ID {}) to expand from queue {}'.format(i, smiles, _id, j))
                    result = self.retroTransformer.get_outcomes(smiles, self.mincount, (self.precursor_prioritization,
                                                                                        self.template_prioritization),
                                                                template_count=self.template_count, 
                                                                mode = self.precursor_score_mode)
                    precursors = result.return_top(n=self.max_branching)
                    self.results_queue.put((_id, smiles, precursors))
                    #print('Worker {} added children of {} (ID {}) to results queue'.format(i, smiles, _id))

                except VanillaQueue.Empty:
                    #print('Queue {} empty for worker {}'.format(j, i))
                    pass
                except Exception as e:
                    print e
            time.sleep(0.01)
            self.idle[i] = True

    def coordinate(self):
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
            except Exception:
                elapsed_time = time.time() - start_time
            
        self.stop()

    def build_tree(self, target):
        self.running = True
        with allow_join_result():
            self.tree_dict[1] = {
                'smiles': target,
                'prod_of': [],
                'rct_of': [],
                'depth': 0,
                'ppg': self.pricer.lookup_smiles(target),
            }
            self.chem_to_id[target] = 1

            if self.max_depth == 1:
                result = self.retroTransformer.get_outcomes(target, self.mincount, (self.precursor_prioritization,
                                                                                    self.template_prioritization),
                                                            template_count=self.template_count,
                                                            mode = self.precursor_score_mode)
                precursors = result.return_top(n=self.max_branching)
                children = self.get_children(precursors)
                self.add_children(children, target, 1)
            else:
                self.prepare()
                self.set_initial_target(target)
                self.coordinate()

    def tree_status(self):
        '''Summarize size of tree after expansion'''
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

    def get_buyable_paths(self, target, max_depth=3, max_branching=25, expansion_time=240, template_prioritization=None,
                          precursor_prioritization=None, nproc=1, mincount=25, chiral=True, max_trees=25, max_ppg=1e10,
                          known_bad_reactions=[], forbidden_molecules=[], template_count=100, precursor_score_mode=gc.max):
        '''Get viable synthesis trees using an iterative deepening depth-first search'''

        self.mincount = mincount
        self.max_depth = max_depth
        self.max_branching = max_branching
        self.expansion_time = expansion_time
        self.template_prioritization = template_prioritization
        self.precursor_prioritization = precursor_prioritization
        self.precursor_score_mode = precursor_score_mode
        self.nproc = nproc
        self.template_count = template_count
        # Load new prices based op specified max price-per-gram
        self.pricer.load(max_ppg=max_ppg)
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

        if not self.celery:
            for i in range(nproc):
                self.idle.append(True)
            if self.max_depth != 1:
                self.expansion_queues = [Queue()
                                         for i in range(self.max_depth - 1)]
            else:
                self.expansion_queues = [Queue()]

        # generate trees
        self.build_tree(target)
        def IDDFS():
            for depth in range(self.max_depth):
                for path in DLS_chem(1, depth, headNode=True):
                    yield chem_dict(1, self.tree_dict[1]['smiles'], self.tree_dict[1]['ppg'], children=path)

        def DLS_chem(chem_id, depth, headNode=False):
            '''Expand at a fixed depth for the current node chem_id.

            headNode indicates whether this is the first (head) node, in which case
            we should expand it even if it is itself buyable'''
            #Copy list so each new branch has separate list.
            if depth == 0:
                # Not allowing deeper - is this buyable?
                if self.tree_dict[chem_id]['ppg']:
                    yield []  # viable node, calling function doesn't need children
            else:
                # Do we need to go deeper?
                if self.tree_dict[chem_id]['ppg'] and not headNode:
                    yield []  # Nope, this is a viable node
                else:
                    # Try going deeper via DLS_rxn function
                    for rxn_id in self.tree_dict[chem_id]['prod_of']:
                        rxn_info_string = ''
                        for path in DLS_rxn(rxn_id, depth):
                            yield [rxn_dict(rxn_id, rxn_info_string, necessary_reagent=self.tree_dict[rxn_id]['necessary_reagent'],
                                            num_examples=self.tree_dict[rxn_id][
                                                'num_examples'], children=path,
                                            template_score=float(self.tree_dict[rxn_id][
                                                                 'template_score']),
                                            smiles='.'.join(sorted([self.tree_dict[x]['smiles'] for x in self.tree_dict[
                                                            rxn_id]['rcts']])) + '>>' + self.tree_dict[chem_id]['smiles'],
                                            tforms=self.tree_dict[rxn_id]['tforms'])]

        def DLS_rxn(rxn_id, depth):
            '''Return children paths starting from a specific rxn_id

            Weird case handling based on 1, 2, or 3 reactants'''

            # Define generators for each reactant node - decrement the depth!
            generators = [DLS_chem(chem_id, depth - 1)
                          for chem_id in self.tree_dict[rxn_id]['rcts']]

            # Only one reactant? easy!
            if len(generators) == 1:
                chem_id = self.tree_dict[rxn_id]['rcts'][0]
                for path in generators[0]:
                    yield [
                        chem_dict(chem_id, self.tree_dict[chem_id][
                                  'smiles'], self.tree_dict[chem_id]['ppg'], children=path)
                    ]

            # Two reactants? want to capture all combinations of each node's
            # options
            elif len(generators) == 2:
                chem_id0 = self.tree_dict[rxn_id]['rcts'][0]
                chem_id1 = self.tree_dict[rxn_id]['rcts'][1]
                for path0 in generators[0]:
                    for path1 in generators[1]:
                        yield [
                            chem_dict(chem_id0, self.tree_dict[chem_id0][
                                      'smiles'], self.tree_dict[chem_id0]['ppg'], children=path0),
                            chem_dict(chem_id1, self.tree_dict[chem_id1][
                                      'smiles'], self.tree_dict[chem_id1]['ppg'], children=path1)
                        ]

            # Three reactants? This is not elegant...
            elif len(generators) == 3:
                chem_id0 = self.tree_dict[rxn_id]['rcts'][0]
                chem_id1 = self.tree_dict[rxn_id]['rcts'][1]
                chem_id2 = self.tree_dict[rxn_id]['rcts'][2]
                for path0 in generators[0]:
                    for path1 in generators[1]:
                        for path2 in generators[2]:
                            yield [
                                chem_dict(chem_id0, self.tree_dict[chem_id0][
                                          'smiles'], self.tree_dict[chem_id0]['ppg'], children=path0),
                                chem_dict(chem_id1, self.tree_dict[chem_id1][
                                          'smiles'], self.tree_dict[chem_id1]['ppg'], children=path1),
                                chem_dict(chem_id2, self.tree_dict[chem_id2][
                                          'smiles'], self.tree_dict[chem_id2]['ppg'], children=path2),
                            ]

            else:
                print('Too many reactants! Only have cases 1-3 programmed')
                raise ValueError(
                    'Too many reactants! Only have cases 1-3 programmed')

        # Generate paths

        import hashlib
        import json
        done_trees = set()
        trees = []
        counter = 0
        for tree in IDDFS():
            hashkey = hashlib.sha1(json.dumps(
                tree, sort_keys=True)).hexdigest()

            # print(tree)
            # print(hashkey)

            if hashkey in done_trees:
                #print('Found duplicate!')
                continue

            done_trees.add(hashkey)
            trees.append(tree)
            counter += 1

            if counter == max_trees:
                MyLogger.print_and_log('Generated {} trees (max_trees met), stopped looking for more...'.format(
                    max_trees), treebuilder_loc)
                break

        return (self.tree_status(), trees)

if __name__ == '__main__':

    MyLogger.initialize_logFile()
    from pymongo import MongoClient
    db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                            'id'], connect=gc.MONGO['connect'])
    TEMPLATE_DB = db_client[gc.RETRO_TRANSFORMS_CHIRAL['database']][
        gc.RETRO_TRANSFORMS_CHIRAL['collection']]
    celery = False
    treedict = []

    treeBuilder = TreeBuilder(celery=celery, mincount=250, mincount_chiral=100)
    # treeBuilder.build_tree('c1ccccc1C(=O)OCCN')
    print treeBuilder.get_buyable_paths('OC(Cn1cncn1)(Cn2cncn2)c3ccc(F)cc3F', max_depth=4, template_prioritization=gc.popularity,
                                        precursor_prioritization=gc.scscore, nproc=16, expansion_time=60, max_trees=500, max_ppg=10,
                                        max_branching = 25,precursor_score_mode =gc.mean)[0]
    print 'done'