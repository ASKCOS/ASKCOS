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

CORRESPONDING_QUEUE = 'tb_worker'

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
def get_top_precursors(smiles, mincount=0, max_branching=20):
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
    precursors = result.return_top(n=max_branching)
    return (smiles, [({'necessary_reagent': precursor['necessary_reagent'],
              'num_examples': precursor['num_examples'],
               'score': precursor['score']}, 
               precursor['smiles_split']) for precursor in precursors])