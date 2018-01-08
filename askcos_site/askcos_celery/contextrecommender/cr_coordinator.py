from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
from celery.result import allow_join_result 
from makeit.synthetic.context.nn_context_recommender import NNContextRecommender
from askcos_site.askcos_celery.contextrecommender.cr_nn_worker import get_n_conditions
import global_config as gc
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
def get_context_recommendations(rxn, n=10, singleSlvt=True, with_smiles=True, context_recommender = ''):
    '''Retrieve a context recommendation given the reaction to attempt.

    rxn = [reacants, products], where each is a list of SMILES
    n = number of contexts to return'''

    global NN_PREDICTOR
    
    print('Context context_recommender worker got a request for rxn {} and n {}'.format(
        rxn, n))
    with allow_join_result():
        if context_recommender == gc.nearest_neighbor:
            res = get_n_conditions.delay(rxn, n=10, singleSlvt=True, with_smiles=True)
            return res.get()
        else:
            raise NotImplementedError