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
from .tb_c_worker import TemplateRelevanceAPIModel, FastFilterAPIModel
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
CORRESPONDING_QUEUE = 'tb_c_worker_preload'
retroTransformer = None


@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    """Configures worker and instantiates RetroTransformer.

    Args:
        options (dict, optional): Used ensure correct queue. (default: {{}})
        **kwargs: Unused.
    """
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A TREE BUILDER WORKER ###')

    # Instantiate and load retro transformer
    global retroTransformer
    retroTransformer = RetroTransformer(template_prioritizer=None, fast_filter=None, use_db=False, load_all=True)
    retroTransformer.load()
    print('### TREE BUILDER WORKER STARTED UP ###')


@shared_task
def get_top_precursors(
        smiles, precursor_prioritizer=None,
        template_set='reaxys', template_prioritizer='reaxys',
        fast_filter=None, max_num_templates=1000,
        max_cum_prob=1, fast_filter_threshold=0.75,
        cluster=True, cluster_method='kmeans', cluster_feature='original',
        cluster_fp_type='morgan', cluster_fp_length=512, cluster_fp_radius=1
    ):
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
        cluster (bool, optional): Whether to cluster results. (default: {True}). This is passed along to RetroResult.return_top()
        cluster_method (str, optional): Clustering method to use ['kmeans', 'hdbscan']. (default: {'kmeans'})
        cluster_feature (str, optional): Features to use for clustering ['original', 'outcomes', 'all']. 'Original' means features that disappear from original target. 'Outcomes' means new features that appear in predicted precursor outcomes. 'All' means the logical 'or' of both. (default: {'original'})
        cluster_fp_type (str, optional): Type of fingerprint to use. Curretnly only 'morgan' is supported. (default: {'morgan'})
        cluster_fp_length (int, optional): Fixed-length folding to use for fingerprint generation. (default: {512})
        cluster_fp_radius (int, optional): Radius to use for fingerprint generation. (default: {1})

    Returns:
        2-tuple of (str, list of dict): SMILES string of input and top
            precursors found.
    """

    template_relevance_hostname = 'template-relevance-{}'.format(template_prioritizer)
    template_prioritizer = TemplateRelevanceAPIModel(
        hostname=template_relevance_hostname, model_name='template_relevance'
    )

    fast_filter_hostname = 'fast-filter'
    fast_filter = FastFilterAPIModel(fast_filter_hostname, 'fast_filter').predict

    global retroTransformer
    result = retroTransformer.get_outcomes(
        smiles, template_set=template_set,
        max_num_templates=max_num_templates, max_cum_prob=max_cum_prob, 
        fast_filter_threshold=fast_filter_threshold, template_prioritizer=template_prioritizer,
        precursor_prioritizer=precursor_prioritizer, fast_filter=fast_filter
    )
    
    return (smiles, result)

@shared_task
def template_relevance(smiles, max_num_templates, max_cum_prob, relevance_model='reaxys'):
    global retroTransformer
    hostname = 'template-relevance-{}'.format(relevance_model)
    template_prioritizer = TemplateRelevanceAPIModel(
        hostname=hostname, model_name='template_relevance'
    )
    scores, indices = template_prioritizer.predict(
        smiles, max_num_templates=max_num_templates, max_cum_prob=max_cum_prob
    )
    if not isinstance(scores, list):
        scores = scores.tolist()
    if not isinstance(indices, list):
        indices = indices.tolist()
    return scores, indices

@shared_task
def apply_one_template_by_idx(*args, **kwargs):
    """Wrapper function for ``RetroTransformer.apply_one_template_by_idx``.

    Returns:
        list of 5-tuples of (int, str, int, list, float): Result of
            applying given template to the molecule.
    """
    global retroTransformer

    template_prioritizer = kwargs.pop('template_prioritizer', 'reaxys')

    hostname = 'template-relevance-{}'.format(template_prioritizer)
    template_prioritizer = TemplateRelevanceAPIModel(
        hostname=hostname, model_name='template_relevance'
    )

    fast_filter_hostname = 'fast-filter'
    fast_filter = FastFilterAPIModel(fast_filter_hostname, 'fast_filter').predict

    kwargs.update({
        'template_prioritizer': template_prioritizer,
        'fast_filter': fast_filter
    })

    return retroTransformer.apply_one_template_by_idx(*args, **kwargs)

@shared_task
def fast_filter_check(*args, **kwargs):
    """Wrapper for fast filter check.

    These workers will already have it initialized. Best way to allow
    independent queries.

    Returns:
        list: Reaction outcomes.
    """
    print('got request for fast filter')
    fast_filter_hostname = 'fast-filter'
    fast_filter = FastFilterAPIModel(fast_filter_hostname, 'fast_filter')
    return fast_filter.predict(*args, **kwargs)
