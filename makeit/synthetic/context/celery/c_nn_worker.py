'''
The role of a context_worker is to take in an attempted reaction and return
a set of conditions to (try to) run the reaction in. These workers rely on 
methods other than the nearest neighbor model to propose contexts.
'''

from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
from synthetic.context.nn_context_recommender import NNContextRecommender
import utilities.i_o.model_loader as model_loader
NN_recommender = None 
CORRESPONDING_QUEUE = 'context_worker'

@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    if 'queues' not in options: 
        return 
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A CONTEXT RECOMMENDER WORKER ###')

    global NN_recommender

    # Get Django settings
    from django.conf import settings

    # Database
    from database import db_client
    db = db_client[settings.INSTANCES['database']]
    INSTANCE_DB = db[settings.INSTANCES['collection']]

    # Setting logging low
    from rdkit import RDLogger
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    
    NN_recommender = model_loader.load_Context_Recommender(max_total_contexts)
    

@shared_task
def get_context_recommendation(rxn, n=1):
    '''
    Retrieve a context recommendation given the reaction to attempt.

    rxn = SMILES
    n = number of contexts to return
    
    '''

    global NN_recommender

    print('Context recommender worker got a request for rxn {} and n {}'.format(
        rxn, n))
   
    return NN_recommender.get_n_conditions(rxn, n=n)