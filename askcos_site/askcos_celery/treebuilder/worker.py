'''
The role of a treebuilder worker is to take a target compound
and apply all retrosynthetic templates to it. The top results
are returned based on the defined mincount and max_branching. The
heuristic chemical scoring function, defined in the transformer
class, is used for prioritization. Each worker pre-loads a
transformer and grabs templates from the database.
'''

from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'tb_worker'
CORRESPONDING_RESERVABLE_QUEUE = 'tb_worker_reservable'

RetroTransformer = None

@celeryd_init.connect
def configure_worker(options={},**kwargs):
    if 'queues' not in options: 
        return 
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A TREE BUILDER WORKER ###')

    global RetroTransformer
    
    # Get Django settings
    from django.conf import settings

    # Database
    from database import db_client
    db = db_client[settings.INSTANCES['database']]
    RETRO_DB = db[settings.RETRO_TRANSFORMS['collection']]

    # Setting logging low
    from rdkit import RDLogger
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    # Import retro transformer class and load
    import makeit.retro.transformer as transformer 
    RetroTransformer = transformer.Transformer(parallel=False, nb_workers=1)
    mincount_retro = settings.RETRO_TRANSFORMS['mincount']
    RetroTransformer.load(RETRO_DB, mincount=mincount_retro, get_retro=True, get_synth=False)
    print('Tree builder worker loaded {} retro templates'.format(RetroTransformer.num_templates))

@shared_task
def get_top_precursors(smiles, mincount=0, max_branching=20, raw_results=False):
    '''Get the precursors for a chemical defined by its SMILES

    smiles = SMILES of node to expand
    mincount = minimum template popularity
    max_branching = maximum number of precursor sets to return, prioritized
        using heuristic chemical scoring function'''

    print('Treebuilder worker was asked to expand {} (mincount {}, branching {})'.format(
        smiles, mincount, max_branching
    ))

    global RetroTransformer

    result = RetroTransformer.perform_retro(smiles, mincount=mincount)
    if result is None:
        precursors = []
    else:
        precursors = result.return_top(n=max_branching)
    if raw_results:
        for i in range(len(precursors)):
            # Must convert ObjectID to string to json-serialize
            precursors[i]['tforms'] = [str(tform) for tform in precursors[i]['tforms']]
        return precursors
    return (smiles, [({'necessary_reagent': precursor['necessary_reagent'],
              'num_examples': precursor['num_examples'],
               'score': precursor['score']}, 
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