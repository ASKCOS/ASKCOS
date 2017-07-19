'''
The role of a treebuilder coordinator is to take a target compound
and build up the retrosynthetic tree by sending individual chemicals
to workers, which each apply the full set of templates. The coordinator
will keep track of the dictionary and ensure unique IDs in addition
to keeping track of the chemical prices using the Pricer module.

The coordinator, finally, returns a set of buyable trees obtained 
from an IDDFS.
'''

from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
from celery.result import allow_join_result 
# NOTE: allow_join_result is only because the treebuilder worker is separate
from celery.exceptions import Terminated
import time

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'tb_coordinator'
Pricer = None


@celeryd_init.connect
def configure_coordinator(options={},**kwargs):
    if 'queues' not in options: 
        return 
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A TREE BUILDER COORDINATOR ###')

    global Pricer
    global get_top_precursors
    
    # Get Django settings
    from django.conf import settings

    # Database
    from database import db_client
    db = db_client[settings.BUYABLES['database']]
    BUYABLE_DB = db[settings.BUYABLES['collection']]
    db = db_client[settings.CHEMICALS['database']]
    CHEMICAL_DB = db[settings.CHEMICALS['collection']]

    ### Prices
    print('Loading prices...')
    import makeit.retro.pricer as pricer
    Pricer = pricer.Pricer()
    Pricer.load(CHEMICAL_DB, BUYABLE_DB)
    print('Loaded known prices')

    # Import worker functions
    from .worker import get_top_precursors, reserve_worker_pool, unreserve_worker_pool
    print('Imported get_top_precursors function from worker')

    import askcos_site.askcos_celery.chiralretro.worker as crworker

    print('Finished initializing treebuilder coordinator')

def chem_dict(_id, smiles, ppg, children=[]):
    '''Chemical object as expected by website'''
    return {
        'id': _id,
        'is_chemical': True,
        'smiles' : smiles,
        'ppg' : ppg,
        'children': children,
    }

def rxn_dict(_id, info, necessary_reagent='', num_examples=0, 
        children=[], smiles='', tforms=[]):
    '''Reaction object as expected by website'''
    return {
        'id': _id,
        'is_reaction': True,
        'info': info,
        'necessary_reagent': necessary_reagent,
        'num_examples': num_examples,
        'children': children,
        'smiles': smiles,
        'tforms': tforms,
    }

def tree_status(tree_dict):
    num_chemicals = 0
    num_reactions = 0
    at_depth = {}
    for _id, node in tree_dict.iteritems():
        depth = node['depth']
        if depth % 1 == 0:
            num_chemicals += 1
        else:
            num_reactions += 1
        if depth not in at_depth:
            at_depth[depth] = 1
        else:
            at_depth[depth] += 1
    return (num_chemicals, num_reactions, at_depth)

def get_trees_iddfs(tree_dict, max_depth, max_trees=25):
    '''Get viable synthesis trees using an iterative deepening depth-first search'''
    from itertools import product

    def IDDFS():
        for depth in range(max_depth+1):
            print('Getting trees (IDDFS starting d={})'.format(depth))
            for path in DLS_chem(1, depth, headNode=True):
                yield chem_dict(1, tree_dict[1]['smiles'], tree_dict[1]['ppg'], children=path)

    def DLS_chem(chem_id, depth, headNode=False):
        '''Expand at a fixed depth for the current node chem_id.

        headNode indicates whether this is the first (head) node, in which case
        we should expand it even if it is itself buyable'''

        if not headNode and tree_dict[chem_id]['smiles'] == tree_dict[1]['smiles']:
            pass # do not let trees terminate in the head node itself
        else:
            if depth == 0:
                # Not allowing deeper - is this buyable?
                if tree_dict[chem_id]['ppg']:
                    yield []  # viable node, calling function doesn't need children
            else:
                # Do we need to go deeper?
                if tree_dict[chem_id]['ppg'] and not headNode:
                    yield [] # Nope, this is a viable node
                else:
                    # Try going deeper via DLS_rxn function
                    for rxn_id in tree_dict[chem_id]['prod_of']:
                        rxn_info_string = ''
                        for path in DLS_rxn(rxn_id, depth):
                            yield [rxn_dict(rxn_id, rxn_info_string, tree_dict[rxn_id]['necessary_reagent'], 
                                 tree_dict[rxn_id]['num_examples'], 
                                 children=path,
                                 smiles='.'.join(sorted([tree_dict[x]['smiles'] for x in tree_dict[rxn_id]['rcts']])) + '>>' + tree_dict[chem_id]['smiles'],
                                 tforms=tree_dict[rxn_id]['tforms'])]

    def DLS_rxn(rxn_id, depth):
        '''Return children paths starting from a specific rxn_id'''
        
        # Define generators for each reactant node - decrement the depth!
        generators = [DLS_chem(chem_id, depth - 1) for chem_id in tree_dict[rxn_id]['rcts']]

        # Only one reactant? easy!
        if len(generators) == 1:
            chem_id = tree_dict[rxn_id]['rcts'][0]
            for path in generators[0]:
                yield [
                    chem_dict(chem_id, tree_dict[chem_id]['smiles'], tree_dict[chem_id]['ppg'], children = path)
                ]

        # Two reactants? want to capture all combinations of each node's options
        elif len(generators) == 2:
            chem_id0 = tree_dict[rxn_id]['rcts'][0]
            chem_id1 = tree_dict[rxn_id]['rcts'][1]
            for path0 in generators[0]:
                for path1 in generators[1]:
                    yield [
                        chem_dict(chem_id0, tree_dict[chem_id0]['smiles'], tree_dict[chem_id0]['ppg'], children = path0),
                        chem_dict(chem_id1, tree_dict[chem_id1]['smiles'], tree_dict[chem_id1]['ppg'], children = path1)
                    ]
                
        # Three reactants? This is not elegant...
        elif len(generators) == 3:
            chem_id0 = tree_dict[rxn_id]['rcts'][0]
            chem_id1 = tree_dict[rxn_id]['rcts'][1]
            chem_id2 = tree_dict[rxn_id]['rcts'][2]
            for path0 in generators[0]:
                for path1 in generators[1]:
                    for path2 in generators[2]:
                        yield [
                            chem_dict(chem_id0, tree_dict[chem_id0]['smiles'], tree_dict[chem_id0]['ppg'], children = path0),
                            chem_dict(chem_id1, tree_dict[chem_id1]['smiles'], tree_dict[chem_id1]['ppg'], children = path1),
                            chem_dict(chem_id2, tree_dict[chem_id2]['smiles'], tree_dict[chem_id2]['ppg'], children = path2),
                        ]

        else:
            print('Too many reactants! Only have cases 1-3 programmed')
            #raise ValueError('Too many reactants! Only have cases 1-3 programmed')
            # do not yield anything

    # Generate paths

    import hashlib 
    import json
    done_trees = set()
    trees = []
    counter = 0
    for tree in IDDFS():

        hashkey = hashlib.sha1(json.dumps(tree, sort_keys=True)).hexdigest()

        if hashkey in done_trees:
            #print('Found duplicate!')
            continue

        done_trees.add(hashkey)
        trees.append(tree)
        counter += 1

        if counter == max_trees:
            print('Generated {} unique trees, stopping looking for more...'.format(max_trees))
            break

    return trees


@shared_task(bind=True)
def get_buyable_paths(self, smiles, mincount=0, max_branching=20, max_depth=3, 
        max_ppg=1e8, max_time=60, max_trees=25, reporting_freq=5,
        known_bad_reactions=[], return_d1_if_no_trees=False, chiral=False):
    '''Get a set of buyable trees for a target compound.

    mincount = minimum template popularity
    max_branching = maximum number of precursor sets to return, prioritized
        using heuristic chemical scoring function
    max_depth = maximum depth to build the tree out to
    max_ppg = maximum price to consider something buyable
    max_time = time for expansion
    reporting_freq = interval (s) for reporting status
    known_bad_reactions = list of reactant smiles for reactions which we
         already know not to work, so we shouldn't propose them

    This function updates its state to provide information about the current
    status of the expansion. Search "on_message task progress" for examples'''

    print('Treebuilder coordinator was asked to expand {}'.format(smiles))
    print('Treebuilder coordinator: mincount {}, max_depth {}, max_branching {}, max_ppg {}, max_time {}, max_trees {}'.format(
        mincount, max_depth, max_branching, max_ppg, max_time, max_trees
    ))

    global Pricer
    from .worker import get_top_precursors, reserve_worker_pool, unreserve_worker_pool
    import askcos_site.askcos_celery.chiralretro.worker as crworker

    tree_dict = {}
    chem_to_id = {smiles: 1}
    buyable_leaves = set()
    pending_results = []

    with allow_join_result(): # required to use .get() from tb_workers

        # Reserve a pool of celery workers
        from celery.exceptions import TimeoutError
        try:
            print('Searching for private worker pool to reserve')
            if chiral:
                private_worker_queue = crworker.reserve_worker_pool.delay().get(timeout=5)
            else:
                private_worker_queue = reserve_worker_pool.delay().get(timeout=5)
            print('Found one! Will use private queue {}'.format(private_worker_queue))
        except TimeoutError:
            raise TimeoutError('Treebuidler coordinator could not find an open pool of workers to reserve!')
        # Wrap everything else in a try/finally block to make sure we unreserve this pool at the end!
        try:
            def expand(smiles, priority=0):
                # Cannot use 'delay' wrapper if we want to use priority arg, so use apply_async
                # Note that we use our private_worker_queue which has been previously reserved
                if chiral:
                    return crworker.get_top_precursors.apply_async(
                        args=(smiles,), 
                        kwargs={'mincount':mincount, 'max_branching':max_branching}, 
                        priority=int(priority),
                        queue=private_worker_queue,
                    )
                return get_top_precursors.apply_async(
                    args=(smiles,), 
                    kwargs={'mincount':mincount, 'max_branching':max_branching}, 
                    priority=int(priority),
                    queue=private_worker_queue,
                )

            # Initialize
            tree_dict[1] = {
                'smiles': smiles,
                'prod_of': [],
                'rct_of': [],
                'depth': 0,
                'ppg': Pricer.lookup_smiles(smiles),
            }
            current_id = 2

            time_start = time.time()
            time_goal = time_start + max_time
            time_last_report = 0.

            # Add first expansion
            pending_results = [expand(smiles, priority=99)]

            while len(pending_results) > 0:
                time_now = time.time()

                # Look for done results
                is_ready = [i for (i, res) in enumerate(pending_results) if res.ready()]
                
                for i in is_ready:
                    (smiles, children) = pending_results[i].get(timeout=1)
                    pending_results[i].forget()
                    _id = chem_to_id[smiles]

                    # Assign unique number
                    for (rxn, mols) in children:

                        # Check if we know this to be a bad reaction
                        # and leave it out of the tree if so
                        rxn_smiles = '.'.join(sorted(mols)) + '>>' + smiles
                        if rxn_smiles in known_bad_reactions:
                            continue

                        rxn_id = current_id
                        current_id += 1 # this is only okay because there is ONE coordinator process

                        # For the parent molecule, record child reactions
                        tree_dict[_id]['prod_of'].append(rxn_id)
                    
                        # For the reaction, keep track of children IDs
                        chem_ids = []
                        for mol in mols:

                            try:
                                chem_id = chem_to_id[mol] # if this does not fail, known

                                # Just need to record that this is a precursor of that rxn
                                tree_dict[chem_id]['rct_of'].append(rxn_id)

                                # Are we less deep than when this was discovered?
                                prev_depth = tree_dict[chem_id]['depth']
                                if tree_dict[_id]['depth'] + 1 < prev_depth:
                                    tree_dict[chem_id]['depth'] = tree_dict[_id]['depth'] + 1

                                    # Do we need to expand now? This is only if
                                    # (a) was previously at max_depth, and (b) not buyable
                                    if prev_depth == max_depth and not tree_dict[chem_id]['ppg']:
                                        # Also, ignore molecules that are way too small
                                        if len(mol) > 4:
                                            pending_results.append(expand(mol, priority=tree_dict[chem_id]['depth']))

                                    # Weird case: the decrease in depth should
                                    # propagate to children, so maybe we need to
                                    # look into the chemical's children and
                                    # figure out if those should be expanded
                                    # again? Not sure how to handle this

                            except KeyError:
                                # New chemical?
                                chem_id = current_id
                                current_id += 1
                                chem_to_id[mol] = chem_id

                                # Check if buyable
                                ppg = Pricer.lookup_smiles(mol, alreadyCanonical=True)
                                if ppg > max_ppg:
                                    ppg = 0.

                                # Add to tree
                                tree_dict[chem_id] = {
                                    'smiles': mol,
                                    'prod_of': [],
                                    'rct_of': [rxn_id],
                                    'depth': tree_dict[_id]['depth'] + 1,
                                    'ppg': ppg
                                }
                                
                                # Check buyability / expandability
                                if ppg:
                                    buyable_leaves.add(chem_id)
                                elif tree_dict[chem_id]['depth'] < max_depth:
                                    # Add to queue to get expanded
                                    pending_results.append(expand(mol, priority=tree_dict[chem_id]['depth']))

                            # Keep track of chem_ids
                            chem_ids.append(chem_id)

                        # Record by overwriting the whole dict value
                        rxn['rcts'] = chem_ids
                        rxn['prod'] = _id
                        rxn['depth'] = tree_dict[_id]['depth'] + 0.5
                        tree_dict[rxn_id] = rxn

                
                # Update list
                pending_results = [res for (i, res) in enumerate(pending_results) if i not in is_ready]

                # Update status
                if time_now - time_last_report >= reporting_freq:
                    (num_chemicals, num_reactions, at_depth) = tree_status(tree_dict)
                    self.update_state(state='EXPANDING', meta={
                        'time_remaining': time_goal - time_now,
                        'num_chemicals': num_chemicals,
                        'num_reactions': num_reactions,
                        'depth_dict': at_depth,
                    })
                    time_last_report = time_now 
                    print('Treebuilder coordinator: updated state: {} chems and {} rxns'.format(num_chemicals, num_reactions))
                    print('at_depth dict: {}'.format(at_depth))
                    print('Treebuilder coordinator: {:.1f} seconds left'.format(time_goal - time_now))

                # Break when appropriate
                if time_now >= time_goal:
                    break

            # Return trees
            print('Getting tree status')
            this_tree_status = tree_status(tree_dict)
            print('Getting trees')
            these_trees = get_trees_iddfs(tree_dict, max_depth, max_trees)
            print('Got {} trees'.format(these_trees))

            # Did we fail to find anything?
            if len(these_trees) == 0 and return_d1_if_no_trees:
                # Return a list of reaction SMILES to try instead
                d1_reactions = []
                for rxn_id in tree_dict[1]['prod_of']:
                    d1_reactions.append((
                        [tree_dict[x]['smiles'] for x in tree_dict[rxn_id]['rcts']],
                        [tree_dict[1]['smiles']],
                        tree_dict[rxn_id]['necessary_reagent'],
                    ))
                return (this_tree_status, these_trees, d1_reactions)
                    
            return (this_tree_status, these_trees)
        
        finally: # should catch KeyboardInterrupt and Terminated
            print('Purging remaining tasks from {}...'.format(private_worker_queue))
            if pending_results != []:
                ## OPTION 1 - REVOKE TASKS, WHICH GETS SENT TO ALL WORKERS REGARDLESS OF TYPE
                #[res.revoke() for res in pending_results]
                # OPTION 2 - DIRECTLY PURGE THE QUEUE (NOTE: HARDCODED FOR AMQP)
                import celery.bin.amqp 
                from askcos_site.celery import app 
                amqp = celery.bin.amqp.amqp(app = app)
                amqp.run('queue.purge', private_worker_queue)

            print('Trying to release private worker pool...')
            if chiral:
                released = crworker.unreserve_worker_pool.apply_async(queue=private_worker_queue, retry=True).get()
            else:
                released = unreserve_worker_pool.apply_async(queue=private_worker_queue, retry=True).get()
            
            if released: 
                print('Released private worker pool!')
            else:
                print('...a pool of workers was lost')