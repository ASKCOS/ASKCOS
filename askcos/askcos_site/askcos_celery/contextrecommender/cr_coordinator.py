from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
from celery.result import allow_join_result
from makeit.synthetic.context.nearestneighbor import NNContextRecommender
from askcos_site.askcos_celery.contextrecommender.cr_nn_worker import get_n_conditions as neighbor_get_n_conditions
from askcos_site.askcos_celery.contextrecommender.cr_network_worker import get_n_conditions as network_get_n_conditions
import makeit.global_config as gc
CORRESPONDING_QUEUE = 'cr_coordinator'
import time


@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A CONTEXT RECOMMENDER COORDINATOR ###')
    print('### CONTEXT RECOMMENDER COORDINATOR STARTED UP ###')


@shared_task
def get_context_recommendations(*args, **kwargs):
    '''Retrieve a context recommendation given the reaction to attempt.

    rxn = [reacants, products], where each is a list of SMILES
    n = number of contexts to return'''

    context_recommender = kwargs.pop('context_recommender', gc.nearest_neighbor)

    print('Context context_recommender worker got a request: {}, {}'.format(args, kwargs))
    with allow_join_result():
        if context_recommender == gc.nearest_neighbor:
            res = neighbor_get_n_conditions.delay(*args, **kwargs)
            return res.get()
        elif context_recommender == gc.neural_network:
            res = network_get_n_conditions.delay(*args, **kwargs)
            return res.get()
        else:
            raise NotImplementedError
        
@shared_task
def get_recommender_types():
    return [gc.nearest_neighbor,gc.neural_network]