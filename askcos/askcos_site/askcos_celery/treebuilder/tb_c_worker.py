"""A worker to build a retrosynthetic tree.

The role of a treebuilder worker is to take a target compound
and apply all retrosynthetic templates to it. The top results
are returned based on the defined mincount and max_branching. The
heuristic chemical scoring function, defined in the transformer
class, is used for prioritization. Each worker pre-loads a
transformer and grabs templates from the database.
"""

from __future__ import absolute_import, unicode_literals, print_function
from django.conf import settings
from celery import shared_task
from celery.signals import celeryd_init
from pymongo import MongoClient
import makeit.global_config as gc
from makeit.retrosynthetic.transformer import RetroTransformer
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
CORRESPONDING_QUEUE = 'tb_c_worker'
CORRESPONDING_RESERVABLE_QUEUE = 'tb_c_worker_reservable'
retroTransformer = None


@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    """Configures worker and instantiates RetroTransformer.

    Args:
        options (dict, optional): Used ensure correct queue?? (default: {{}})
        **kwargs: Unused.
    """
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A TREE BUILDER WORKER ###')

    # Instantiate and load retro transformer
    global retroTransformer
    retroTransformer = RetroTransformer(celery=True)

    retroTransformer.load(chiral=True)
    print(retroTransformer.fast_filter.evaluate('CCCCCCO.CCCCBr', 'CCCCCCOCCCC'))
    print('### TREE BUILDER WORKER STARTED UP ###')


@shared_task
def get_top_precursors(smiles, template_prioritizer, precursor_prioritizer, mincount=0,
                       max_branching=20, template_count=10000, mode=gc.max, max_cum_prob=1, apply_fast_filter=False, filter_threshold=0.8):
    """Get the precursors for a chemical defined by its SMILES.

    Args:
        smiles (str): SMILES of node to expand.
        template_prioritizer (str): Keyword for which prioritization method for
            the templates should be used. Keywords can be found in
            global_config.
        precursor_prioritizer (str): Keyword for which prioritization method for
            the precursors should be used.
        mincount (int, optional): Minimum template popularity. (default: {0})
        max_branching (int, optional): Maximum number of precursor sets to
            return, prioritized using heuristic chemical scoring function.
            (default: {20})
        template_count (int, optional): Maximum number of templates to consider.
            (default: {10000})
        mode (str, optional): Mode to use for merging list of scores. Used by
            prioritizers. (default: {gc.max})
        max_cum_prob (float, optional): Maximum cumulative probability of
            selected relevant templates. (default: {1})
        apply_fast_filter (bool, optional): Whether to use the fast filter to
            filter precursors. (default: {False})
        filter_threshold (float, optional): Threshold to use for fast filter.
            (default: {0.8})

    Returns:
        2-tuple of (str, list of dict): SMILES string of input and top
            precursors found.
    """

    # print('Treebuilder worker was asked to expand {} (mincount {}, branching {}) using {} and {}'.format(
    #    smiles, mincount, max_branching, template_prioritizer, precursor_prioritizer
    #))

    global retroTransformer
    result = retroTransformer.get_outcomes(
        smiles, mincount, (precursor_prioritizer,
                           template_prioritizer), template_count=template_count, mode=mode,
        max_cum_prob=max_cum_prob, apply_fast_filter=apply_fast_filter, filter_threshold=filter_threshold)

    # print(result)

    precursors = result.return_top(n=max_branching)
    return (smiles, precursors)

@shared_task
def apply_one_template_by_idx(*args, **kwargs):
    """Wrapper function for RetroTransformer.apply_one_template_by_idx."""
    return retroTransformer.apply_one_template_by_idx(*args, **kwargs)

@shared_task
def fast_filter_check(*args, **kwargs):
    """Wrapper for fast filter check.

    These workers will already have it initialized. Best way to allow
    independent queries.
    """
    print('got request for fast filter')
    return retroTransformer.fast_filter.evaluate(*args, **kwargs)


@shared_task(bind=True)
def reserve_worker_pool(self):
    """Reserves pool of workers.

    Called by a tb_coordinator to reserve this pool of workers to do a tree
    expansion. This is accomplished by changing what queue(s) this pool
    listens to.
    """
    hostname = self.request.hostname
    private_queue = CORRESPONDING_QUEUE + '_' + hostname
    print('Tried to reserve this worker!')
    print('I am {}'.format(hostname))
    print('Telling myself to ignore the {} and {} queues'.format(
        CORRESPONDING_QUEUE, CORRESPONDING_RESERVABLE_QUEUE))
    from askcos_site.celery import app
    app.control.cancel_consumer(CORRESPONDING_QUEUE, destination=[hostname])
    app.control.cancel_consumer(
        CORRESPONDING_RESERVABLE_QUEUE, destination=[hostname])

    # *** purge the queue in case old jobs remain
    import celery.bin.amqp
    amqp = celery.bin.amqp.amqp(app=app)
    amqp.run('queue.purge', private_queue)
    print('Telling myself to only listen to the new {} queue'.format(private_queue))
    app.control.add_consumer(private_queue, destination=[hostname])
    return private_queue


@shared_task(bind=True)
def unreserve_worker_pool(self):
    """Releases this worker pool so it can listen to the original queues."""
    hostname = self.request.hostname
    private_queue = CORRESPONDING_QUEUE + '_' + hostname
    print('Tried to unreserve this worker!')
    print('I am {}'.format(hostname))
    print('Telling myself to ignore the {} queue'.format(private_queue))
    from askcos_site.celery import app
    app.control.cancel_consumer(private_queue, destination=[hostname])
    print('Telling myself to only listen to the {} and {} queues'.format(
        CORRESPONDING_QUEUE, CORRESPONDING_RESERVABLE_QUEUE))
    app.control.add_consumer(CORRESPONDING_QUEUE, destination=[hostname])
    app.control.add_consumer(
        CORRESPONDING_RESERVABLE_QUEUE, destination=[hostname])
    return True
