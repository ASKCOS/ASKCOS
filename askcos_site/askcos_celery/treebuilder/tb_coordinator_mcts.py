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
from collections import defaultdict
from celery.result import allow_join_result
# NOTE: allow_join_result is only because the treebuilder worker is separate
from celery.exceptions import Terminated
import time
from rdkit import RDLogger
import makeit.global_config as gc
from makeit.utilities.buyable.pricer import Pricer
from makeit.retrosynthetic.mcts.tree_builder import MCTS as MCTSTreeBuilder
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'tb_coordinator_mcts'


@celeryd_init.connect
def configure_coordinator(options={}, **kwargs):
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A TREE BUILDER MCTS COORDINATOR ###')

    global treeBuilder
    global evaluator

    treeBuilder = MCTSTreeBuilder(celery=True, nproc=8) # 8 active pathways
    print('Finished initializing treebuilder MCTS coordinator')


@shared_task()
def get_buyable_paths(*args, **kwargs):
    print('Treebuilder MCTS coordinator was asked to expand {}'.format(args[0]))
    result = treeBuilder.get_buyable_paths(*args, **kwargs)
    print('Task completed, returning results.')
    return result
