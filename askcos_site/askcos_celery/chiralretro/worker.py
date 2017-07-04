'''
The role of a forward predictor worker is to apply a subset of 
templates and generate candidate edits
'''

from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
import itertools
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'cr_worker'
RetroTransformer = None
apply_one_retrotemplate = None
rdchiralReactants = None

@celeryd_init.connect
def configure_worker(options={},**kwargs):
    if 'queues' not in options: 
        return 
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A CHIRAL RETRO WORKER ###')

    global RetroTransformer
    from .common import get_retro_transformer
    RetroTransformer = get_retro_transformer()
    print('Got RetroTransformer')

    global apply_one_retrotemplate
    from makeit.webapp.transformer_v3 import apply_one_retrotemplate
    print('Imported transformer_v3')
    
    global rdchiralReactants
    from rdchiral.main import rdchiralReactants
    print('Finished initializing chiral retro worker')

@shared_task
def get_chiral_precursor_batch(smiles, start_at, end_at):
    '''Apply chiral retro templates to a target smiles string. We
    use chunks (start_at, end_at) to have fewer queue messages

    smiles = SMILES of target
    start_at = index of templates to start at
    end_at = index of templates to end at'''

    global RetroTransformer
    global rdchiralReactants

    rct = rdchiralReactants(smiles)

    # Each template returns a list, sometimes an empty list
    precursors = itertools.chain.from_iterable(
        map(lambda x: apply_one_retrotemplate(rct, smiles, x, return_as_tup=True), 
        RetroTransformer.templates[start_at:end_at])
    )
    # Don't return empty lists where templates did not apply
    return [precursor for precursor in precursors if precursor]