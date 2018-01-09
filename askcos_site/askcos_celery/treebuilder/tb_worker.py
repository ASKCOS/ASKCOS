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
from pymongo import MongoClient
import makeit.global_config as gc 
from makeit.retro_synthetic.retro_transformer import RetroTransformer
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
CORRESPONDING_QUEUE = 'tb_worker'
CORRESPONDING_RESERVABLE_QUEUE = 'tb_worker_reservable'

@celeryd_init.connect
def configure_worker(options={},**kwargs):
    
    if 'queues' not in options: 
        return 
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A TREE BUILDER WORKER ###')
    # Get Django settings
    from django.conf import settings
    
    # Database
    db_client = MongoClient('mongodb://guest:guest@rmg.mit.edu:27017/admin', serverSelectionTimeoutMS = 2000, connect=True)
    db = db_client[settings.RETRO_TRANSFORMS['database']]
    RETRO_DB = db[settings.RETRO_TRANSFORMS['collection']]
    
    global retroTransformer
    
    #Instantiate and load retro transformer
    retroTransformer = RetroTransformer(TEMPLATE_DB = RETRO_DB, mincount = settings.RETRO_TRANSFORMS['mincount'])
    
    retroTransformer.load()
    print('### TREE BUILDER WORKER STARTED UP ###')


@shared_task
def get_top_precursors(smiles, template_prioritizer, precursor_prioritizer, mincount=25, max_branching=20):
    '''Get the precursors for a chemical defined by its SMILES

    smiles = SMILES of node to expand
    mincount = minimum template popularity
    max_branching = maximum number of precursor sets to return, prioritized
        using heuristic chemical scoring function
    chiral = whether or not to use the version of the transformer that takes chriality into account.
    template_prioritizer = keyword for which prioritization method for the templates should be used, keywords can be found in global_config
    precursor_prioritizer = keyword for which prioritization method for the precursors should be used.'''

    print('Treebuilder worker was asked to expand {} (mincount {}, branching {})'.format(
        smiles, mincount, max_branching
    ))
<<<<<<< HEAD:askcos_site/askcos_celery/treebuilder/tb_worker.py
    
    global retroTransformer
    result = retroTransformer.get_outcomes(smiles, mincount, (precursor_prioritizer, template_prioritizer))
    precursors = result.return_top(n=max_branching)
    print('Task completed, returning results.')
    return (smiles, precursors)
=======

    global RetroTransformer

    result = RetroTransformer.perform_retro(smiles, mincount=mincount)
    if result is None:
        precursors = []
    else:
        precursors = result.return_top(n=max_branching)
    
    # Must convert ObjectID to string to json-serialize
    for i in range(len(precursors)):
        precursors[i]['tforms'] = [str(tform) for tform in precursors[i]['tforms']]

    if raw_results:
        return precursors
        
    return (smiles, [({'necessary_reagent': precursor['necessary_reagent'],
              'num_examples': precursor['num_examples'],
               'score': precursor['score'],
               'tforms': precursor['tforms']}, 
               precursor['smiles_split']) for precursor in precursors])
>>>>>>> 53b125d4bd2cee599ca71a780a384abe74751830:askcos_site/askcos_celery/treebuilder/worker.py

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