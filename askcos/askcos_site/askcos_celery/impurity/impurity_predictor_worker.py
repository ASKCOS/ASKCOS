from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger
import time
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

WLN_predictor = None
transformer_predictor = None
CORRESPONDING_QUEUE = 'atom_mapping_worker'

@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    print(options)
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A IMPURITY PREDICTOR WORKER ###')
    from makeit.synthetic.evaluation.template_free import TemplateFreeNeuralNetScorer
    global WLN_predictor
    global transformer_predictor

    try:
        WLN_predictor = TemplateFreeNeuralNetScorer()
        transformer_predictor = None
    except Exception as e:
        print(e)
        raise (e)
    print('Initialized')


@shared_task
def predict_reaction(reactants_smiles, predictor='WLN forward predictor'):
    """
    Args:
        reactants_smiles:
        predictor:
    Returns:

    """
    global WLN_predictor
    global transformer_predictor

    outcome = []
    if predictor == 'WLN forward predictor':
        try:
            outcome = WLN_predictor.evaluate(reactants_smiles)
        except Exception as e:
            print(e)
    if predictor == 'Molecular transformer':
        try:
            outcome = transformer_predictor.evaluate(reactants_smiles)
        except Exception as e:
            print(e)

    return outcome
