'''
The role of a forward predictor coordinator is to generate a plausibility
score for an intended reaction. It offloads the template-expansion
work onto a separate worker pool. Returned candidate edits are used to
build up the tensors used with the trained Keras model for forward
prediction.
'''

from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
from celery.result import allow_join_result 
# NOTE: allow_join_result is only because the worker is separate
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
import rdkit.Chem as Chem

CORRESPONDING_QUEUE = 'cr_coordinator'
template_counts = None
RetroTransformer = None

@celeryd_init.connect
def configure_worker(options={},**kwargs):
    if 'queues' not in options: 
        return 
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A CHIRAL RETRO COORDINATOR ###')

    global template_counts
    global RetroTransformer
    from .common import get_retro_transformer
    RetroTransformer = get_retro_transformer()
    template_counts = [template['count'] for template in RetroTransformer.templates]

    print('Finished initializing chiral retro coordinator')

@shared_task
def get_top_chiral_precursors(smiles, mincount=0, max_branching=20, 
        raw_results=False, chunksize=2000):
    '''Get the precursors for a chemical defined by its SMILES. This
    function is meant for CHIRAL retro steps. Since these are very slow, 
    the tempaltes get parallelized over a number of workers.

    smiles = SMILES of node to expand
    mincount = minimum template popularity
    max_branching = maximum number of precursor sets to return, prioritized
        using heuristic chemical scoring function'''

    print('Chiral retro coordinator was asked to expand {} (mincount {}, branching {})'.format(
        smiles, mincount, max_branching
    ))

    global RetroTransformer

    global template_counts
    from .worker import get_chiral_precursor_batch

    # Make sure we can parse
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    
    # Figure out what templates we are going to use
    end_at = len(template_counts)
    for (i, count) in enumerate(template_counts):
        if count < mincount:
            end_at = i
            break

    # Chunk and add to queue
    pending_results = []
    for start_at in range(0, end_at, chunksize):
        pending_results.append(
            get_chiral_precursor_batch.delay(smiles, start_at, start_at + chunksize)
        )
    print('Added {} chunks'.format(len(pending_results)))

    # Initialize result
    from makeit.webapp.transformer_v3 import RetroResult, RetroPrecursor
    result = RetroResult(smiles)
    add_precursor = result.add_precursor

    # Look for results
    while len(pending_results) > 0:
        is_ready = [i for (i, res) in enumerate(pending_results) if res.ready()]
        with allow_join_result(): # required to use .get()
            for i in is_ready:
                for precursor in pending_results[i].get(timeout=1):
                    # Add precursor (after converting to RetroPrecursor from tuple)
                    add_precursor(RetroPrecursor(
                        smiles_list=precursor[0],
                        template_id=precursor[1],
                        num_examples=precursor[2],
                        necessary_reagent=precursor[3]
                    ), Pricer=RetroTransformer.Pricer)
                pending_results[i].forget()
        # Update list
        pending_results = [res for (i, res) in enumerate(pending_results) if i not in is_ready]

    # Done! return results in whatever format is needed
    precursors = result.return_top(n=max_branching)
    if raw_results:
        #for i in range(len(precursors)):
            ## Must convert ObjectID to string to json-serialize
            #precursors[i]['tforms'] = [str(tform) for tform in precursors[i]['tforms']]
        # ---> Will already be strings since the workers had to convert them
        return precursors
    return (smiles, [({'necessary_reagent': precursor['necessary_reagent'],
              'num_examples': precursor['num_examples'],
               'score': precursor['score'],
               'tforms': precursor['tforms']}, 
               precursor['smiles_split']) for precursor in precursors])