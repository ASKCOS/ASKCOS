from __future__ import absolute_import, unicode_literals, print_function
from django.conf import settings
from celery import shared_task
from celery.signals import celeryd_init
from pymongo import MongoClient
from celery.result import allow_join_result
from celery.exceptions import Terminated
import time
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'sc_coordinator'


@celeryd_init.connect
def configure_coordinator(options={}, **kwargs):
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A SCORING COORDINATOR ###')
    from makeit.synthetic.evaluation.evaluator import Evaluator

    global evaluator

    evaluator = Evaluator(celery=True)
    #evaluator.evaluate('CCC(=O)O.CCCCN','CCC(=O)NCCCC',[(65, 'CCO', '', '', 24, 85)], mincount=100,forward_scorer='templatebased')
    print('### SCORING COORDINATOR STARTED UP ###')


@shared_task(bind=True)
def evaluate(self, reactant_smiles, target, contexts, **kwargs):
    print('Scoring Coordinator was asked to evaluate {} to {}'.format(
        reactant_smiles, target))
    result = evaluator.evaluate(reactant_smiles, target, contexts, **kwargs)
    print('Task completed, returning results.')
    print(result)
    return result
