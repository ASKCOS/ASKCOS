"""
The role of a context_worker is to take in an attempted reaction and return
a set of conditions to (try to) run the reaction in. Each worker will
load a pre-trained neural network model. For each request, this worker
must query the database to get details about the instance.
"""

from __future__ import absolute_import, unicode_literals, print_function
from django.conf import settings
from celery import shared_task
from celery.signals import celeryd_init
import makeit.global_config as gc

CORRESPONDING_QUEUE = 'cr_network_worker'


@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A NEURAL NETWORK CONTEXT RECOMMENDER WORKER ###')
    from makeit.synthetic.context.neuralnetwork import NeuralNetContextRecommender

    global recommender

    # Setting logging low
    from rdkit import RDLogger
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    try:
        recommender = NeuralNetContextRecommender()
        recommender.load()
    except Exception as e:
        print(e)
    print('Loaded context recommendation model')

    print('### NEURAL NETWORK CONTEXT RECOMMENDER STARTED UP ###')


@shared_task
def get_n_conditions(*args, **kwargs):
    """Retrieve a context recommendation given the reaction to attempt.

    rxn = [reacants, products], Where each is a list of SMILES.
    n = Number of contexts to return.
    """

    print('Context recommender worker got a request: {}, {}'.format(args, kwargs))
    res = recommender.get_n_conditions(*args, **kwargs)

    print('Task completed, returning results.')
    return res
