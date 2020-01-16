from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger
import time
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

reaxys_inspector = None
pistachio_inspector = None
CORRESPONDING_QUEUE = 'atom_mapping_worker'

@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    print(options)
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A IMPURITY PREDICTOR WORKER ###')
    from makeit.synthetic.evaluation.fast_filter import FastFilerImpurityInspector
    global reaxys_inspector
    global pistachio_inspector

    try:
        reaxys_inspector = FastFilerImpurityInspector()
        pistachio_inspector = None
    except Exception as e:
        print(e)
        raise (e)
    print('Initialized')


@shared_task
def inspect_reaction(rxnsmiles, inspector='Reaxys inspector'):
    """
    Args:
        rxnsmiles:
        inspector: Reaxys inspector OR Pistachio inspector
    Returns:

    """
    global reaxys_inspector
    global pistachio_inspector

    score = 0
    if inspector == 'Reaxys inspector':
        try:
            score = reaxys_inspector.evaluate(rxnsmiles)
        except Exception as e:
            print(e)
    if inspector == 'Pistachio inspector':
        try:
            score = pistachio_inspector.evaluate(rxnsmiles)
        except Exception as e:
            print(e)

    return score
