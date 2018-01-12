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
from django.conf import settings
from celery import shared_task
from celery.signals import celeryd_init
from pymongo import MongoClient
from celery.result import allow_join_result
# NOTE: allow_join_result is only because the treebuilder worker is separate
from celery.exceptions import Terminated
import time
from askcos_site.askcos_celery.treebuilder.tb_worker import get_top_precursors, reserve_worker_pool, unreserve_worker_pool
from rdkit import RDLogger
from makeit.utilities.buyable.pricer import Pricer
from makeit.retrosynthetic.tree_builder import TreeBuilder
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'tb_coordinator'

from makeit.synthetic.forward_evaluation.evaluator import Evaluator


@celeryd_init.connect
def configure_coordinator(options={}, **kwargs):
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A TREE BUILDER COORDINATOR ###')

    global treeBuilder
    global evaluator

    evaluator = Evaluator(celery=True)

    # Database
    from database import db_client
    db = db_client[settings.BUYABLES['database']]
    BUYABLE_DB = db[settings.BUYABLES['collection']]
    db = db_client[settings.CHEMICALS['database']]
    CHEMICAL_DB = db[settings.CHEMICALS['collection']]

    # Prices
    print('Loading prices...')
    pricer = Pricer(CHEMICALS=CHEMICAL_DB, BUYABLES=BUYABLE_DB)
    print('Loaded known prices')
    treeBuilder = TreeBuilder(celery=True, pricer=pricer)

    print('Finished initializing treebuilder coordinator')


@shared_task(bind=True)
def get_buyable_paths(self, smiles, template_prioritization, precursor_prioritization, mincount=0, max_branching=20,
                      max_depth=3, max_ppg=1e8, max_time=60, max_trees=25, reporting_freq=5, known_bad_reactions=[],
                      return_d1_if_no_trees=False, chiral=False):
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
    result = treeBuilder.get_buyable_paths(smiles, max_depth=max_depth, max_branching=max_branching, expansion_time=max_time,
                                           template_prioritization=template_prioritization, precursor_prioritization=precursor_prioritization,
                                           mincount=mincount, chiral=chiral, max_trees=max_trees, max_ppg=max_ppg, known_bad_reactions=known_bad_reactions)
    print('Task completed, returning results.')
    return result
