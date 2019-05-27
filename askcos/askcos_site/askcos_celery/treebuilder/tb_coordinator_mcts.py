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
import celery
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
from django.db import models
from askcos_site.main.models import SavedResults
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'tb_coordinator_mcts'

client = MongoClient(
    gc.MONGO['path'],
    gc.MONGO['id'],
    connect=gc.MONGO['connect']
)
results_db = client['results']
results_collection = results_db['results']

def update_result_state(id_, state):
    result = SavedResults.objects.get(result_id=id_)
    result.result_state = state
    result.save()
    return

def save_results(result, settings, task_id):
    doc = {
        '_id': task_id,
        'result': result,
        'settings': settings
    }
    results_collection.insert_one(doc)

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


@shared_task(trail=False)
def get_buyable_paths(*args, **kwargs):
    print('Treebuilder MCTS coordinator was asked to expand {}'.format(args[0]))
    _id = get_buyable_paths.request.id
    try:
        result = treeBuilder.get_buyable_paths(*args, **kwargs)
    except:
        update_result_state(_id, 'failed')
        return {}
    update_result_state(_id, 'completed')
    settings = {'smiles': args[0]}
    settings.update(kwargs)
    save_results(result, settings, _id)
    print('Task completed, returning results.')
    return result
