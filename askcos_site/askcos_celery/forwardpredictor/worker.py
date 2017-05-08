'''
The role of a forward predictor worker is to apply a subset of 
templates and generate candidate edits
'''

from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
import time

CORRESPONDING_QUEUE = 'fp_worker'
templates = None

@celeryd_init.connect
def configure_worker(options={},**kwargs):
    if 'queues' not in options: 
        return 
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A FORWARD PREDICTOR WORKER ###')

    global templates
    
    # Get Django settings
    from django.conf import settings

    # Database
    from database import db_client
    db = db_client[settings.SYNTH_TRANSFORMS['database']]
    SYNTH_DB = database[settings.SYNTH_TRANSFORMS['collection']]

    # Load templates
    from .common import load_templates
    mincount_synth = settings.SYNTH_TRANSFORMS['mincount']
    templates = load_templates(SYNTH_DB=SYNTH_DB, mincount=mincount_synth)
    print('Finished initializing forward predictor worker')

@shared_task
def get_candidate_edits(reactants, start_at, end_at):
    '''Evaluate the plausibility of a proposed forward reaction

    reactants = SMILES of reactants
    start_at = index of templates to start at
    end_at = index of templates to end at'''

    print('Forward predictor worker was asked to expand {} ({}->{})'.format(reactants, start_at, end_at))

    global templates

    candidate_edits = []
    for i in range(start_at, end_at):
        try:


        except IndexError: # out of range!
            break
    return candidate_edits