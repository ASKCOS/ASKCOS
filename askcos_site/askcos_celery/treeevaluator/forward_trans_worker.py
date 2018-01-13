'''
The role of a forward predictor worker is to apply a subset of 
templates and generate candidate edits
'''

from __future__ import absolute_import, unicode_literals, print_function
from django.conf import settings
from celery import shared_task
from celery.signals import celeryd_init
import time
import rdkit.Chem as Chem
from makeit.synthetic.enumeration.transformer import ForwardTransformer
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'ft_worker'
templates = None


@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    if 'queues' not in options:
        print('Queues not in options')
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        print('Queues not in options as key')
        return
    print('### STARTING UP A FORWARD ENUMERATION WORKER ###')

    global templates
    global forwardTransformer

    # Database
    from database import db_client
    db = db_client[settings.SYNTH_TRANSFORMS['database']]
    SYNTH_DB = db[settings.SYNTH_TRANSFORMS['collection']]
    # Load templates
    mincount_synth = settings.SYNTH_TRANSFORMS['mincount']
    forwardTransformer = ForwardTransformer(
        celery=True, TEMPLATE_DB=SYNTH_DB, mincount=mincount_synth)
    forwardTransformer.load()
    print('### FORWARD ENUMERATION WORKER STARTED UP ###')


@shared_task
def get_outcomes(reactants_smiles, mincount, start_at, end_at, template_prioritization, template_count=10000):
    '''Apply forward templates to a atom-mapped reactant pool. We
    use chunks (start_at, end_at) to have fewer queue messages

    reactants_smiles = SMILES of reactants (atom-mapped)
    start_at = index of templates to start at
    end_at = index of templates to end at'''
    print('Forward transformer worker was asked to expand {} (mincount {}, from {} to {})'.format(
        reactants_smiles, mincount, start_at, end_at))
    global forwardTransformer
    print('Task completed, returning results.')
    return forwardTransformer.get_outcomes(reactants_smiles, mincount, template_prioritization,
                                           start_at=start_at, end_at=end_at, template_count=template_count)


@shared_task
def template_count():
    global forwardTransformer
    return forwardTransformer.template_count()
