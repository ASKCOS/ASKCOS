from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger
import time
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

wln_mapper = None
heuristic_mapper = None
CORRESPONDING_QUEUE = 'atom_mapping_worker'

@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    print(options)
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A IMPURITY PREDICTOR WORKER ###')
    from makeit.synthetic.atom_mapper.wln_mapper import WLN_AtomMapper
    global wln_mapper
    global heuristic_mapper

    try:
        wln_mapper = WLN_AtomMapper()
        heuristic_mapper = None
    except Exception as e:
        print(e)
        raise (e)
    print('Initialized')


@shared_task
def get_atom_mapping(rxnsmiles, mapper='WLN atom mapper'):
    """
    Args:
        rxnsmiles:
        mapper: two options: WLN atom mapper & Heuristic mapper
    Returns:

    """
    global wln_mapper
    global heuristic_mapper

    rxnsmiles_mapped = ''
    if mapper == 'WLN atom mapper':
        try:
            rxnsmiles_mapped = wln_mapper.evaluate(rxnsmiles)
        except Exception as e:
            print(e)
    if mapper == 'Heuristic mapper':
        try:
            rxnsmiles_mapped = heuristic_mapper.evaluate(rxnsmiles)
        except Exception as e:
            print(e)

    if not rxnsmiles_mapped:
        print('Failed to map the given reaction smiles')

    return rxnsmiles_mapped
