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
CORRESPONDING_RESERVABLE_QUEUE = 'cr_worker_reservable'


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


@shared_task
def get_top_precursors(smiles, mincount=0, max_branching=20, raw_results=False):
    '''Get the precursors for a chemical defined by its SMILES

    smiles = SMILES of node to expand
    mincount = minimum template popularity
    max_branching = maximum number of precursor sets to return, prioritized
        using heuristic chemical scoring function'''

    print('Chiral retro worker was asked to expand {} (mincount {}, branching {})'.format(
        smiles, mincount, max_branching
    ))

    global RetroTransformer

    result = RetroTransformer.perform_retro(smiles, mincount=mincount)
    if result is None:
        precursors = []
    else:
        precursors = result.return_top(n=max_branching)

    for i in range(len(precursors)):
        precursors[i]['tforms'] = [str(tform) for tform in precursors[i]['tforms']]

    if raw_results:
        return precursors
    return (smiles, [({'necessary_reagent': precursor['necessary_reagent'],
              'num_examples': precursor['num_examples'],
               'score': precursor['score'],
               'tforms': precursor['tforms']}, 
               precursor['smiles_split']) for precursor in precursors])

@shared_task(bind=True)
def reserve_worker_pool(self):
    '''Called by a tb_coordinator to reserve this
    pool of workers to do a tree expansion. This is
    accomplished by changing what queue(s) this pool
    listens to'''
    hostname = self.request.hostname
    private_queue = CORRESPONDING_QUEUE + '_' + hostname
    print('Tried to reserve this worker!')
    print('I am {}'.format(hostname))
    print('Telling myself to ignore the {} and {} queues'.format(CORRESPONDING_QUEUE, CORRESPONDING_RESERVABLE_QUEUE))
    from askcos_site.celery import app 
    app.control.cancel_consumer(CORRESPONDING_QUEUE, destination=[hostname])
    app.control.cancel_consumer(CORRESPONDING_RESERVABLE_QUEUE, destination=[hostname])

    # *** purge the queue in case old jobs remain
    import celery.bin.amqp 
    amqp = celery.bin.amqp.amqp(app = app)
    amqp.run('queue.purge', private_queue)
    print('Telling myself to only listen to the new {} queue'.format(private_queue))
    app.control.add_consumer(private_queue, destination=[hostname])
    return private_queue

@shared_task(bind=True)
def unreserve_worker_pool(self):
    '''Releases this worker pool so it can listen
    to the original queues'''
    hostname = self.request.hostname
    private_queue = CORRESPONDING_QUEUE + '_' + hostname
    print('Tried to unreserve this worker!')
    print('I am {}'.format(hostname))
    print('Telling myself to ignore the {} queue'.format(private_queue))
    from askcos_site.celery import app 
    app.control.cancel_consumer(private_queue, destination=[hostname])
    print('Telling myself to only listen to the {} and {} queues'.format(CORRESPONDING_QUEUE, CORRESPONDING_RESERVABLE_QUEUE))
    app.control.add_consumer(CORRESPONDING_QUEUE, destination=[hostname])
    app.control.add_consumer(CORRESPONDING_RESERVABLE_QUEUE, destination=[hostname])
    return True