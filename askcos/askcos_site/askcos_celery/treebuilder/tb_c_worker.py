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
from makeit.utilities.fingerprinting import create_rxn_Morgan2FP_separately
from rdkit import RDLogger, Chem
from rdkit.Chem import AllChem
import numpy as np
from bson import ObjectId
from scipy.special import softmax
import rdkit.Chem as Chem
from rdkit.Chem import AllChem

from ..tfserving import TFServingAPIModel

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
CORRESPONDING_QUEUE = 'tb_c_worker'
CORRESPONDING_RESERVABLE_QUEUE = 'tb_c_worker_reservable'
retroTransformer = None


class TemplateRelevanceAPIModel(TFServingAPIModel):
    """Template relevance Tensorflow API Model. Overrides input and output transformation methods with template relevance specific methods.

    Attributes:
        hostname (str): hostname of service serving tf model.
        model_name (str): Name of model provided to tf serving.
        version (int): version of the model to use when serving
    """
    def transform_input(self, smiles, fp_length=2048, fp_radius=2, **kwargs):
        """Transforms the input for the API model from a SMILES string to a fingerprint bit vector

        Args:
            smiles (str): SMILES string of input molecule
            fp_length (int): Length of desired fingerprint. Must agree with model input shape
            fp_radius (int): Radius of desired fingerprint. Should agree with parameter sed when model was trained

        Returns:
            list: Fingerprint bit vector as a list
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return np.zeros((self.fp_length,), dtype=np.float32)
        return np.array(
            AllChem.GetMorganFingerprintAsBitVect(
                mol, fp_radius, nBits=fp_length, useChirality=True
            ), dtype=np.float32
        ).reshape(1, -1).tolist()
    
    def transform_output(self, pred, max_num_templates=100, max_cum_prob=0.995, **kwargs):
        """Transforms output of API model to return the top scores and indices for the output classes (templates)

        Args:
            pred (np.array): numpy array of output from prediction
            max_num_templates (int): Maximum number of template scores/indices to return from the prediction
            max_cum_prob (float): Maximum cumulative probability of templates to be returned. This offers another convenient way to limit the number of templates returned to those are have been given high scores
        """
        scores = softmax(pred)
        indices = np.argsort(-scores)[:max_num_templates]
        scores = scores[indices]
        cum_scores = np.cumsum(scores)
        if max_cum_prob >= cum_scores[-1]:
            truncate = -1
        else:
            truncate = np.argmax(cum_scores > max_cum_prob)
        return scores[:truncate], indices[:truncate]


class FastFilterAPIModel(TFServingAPIModel):
    """Fast filter Tensorflow API Model. Overrides input and output transformation methods with fast filter specific methods.

    Attributes:
        hostname (str): hostname of service serving tf model.
        model_name (str): Name of model provided to tf serving.
        version (int): version of the model to use when serving
    """
    def transform_input(self, reactant_smiles, target, rxnfpsize=2048, pfpsize=2048, useFeatures=False):
        """Transforms the input for the API model from SMILES strings to product and reaction fingerprints

        Args:
            reactant_smiles (str): SMILES string of the reactants
            target (str): SMILES string of the target
            rxnfpsize (int): Length of desired fingerprint for the reaction. Must agree with model input shape
            pfpsize (int): Length of desired fingerprint for the product. Must agree with model input shape
            useFeatures (bool): Flag to use features or not when generating fingerprint. Should agree with how model was trained

        Returns:
            list of dict: Input fingerprints, formatted for a call to the tensorflow API model
        """
        pfp, rfp = create_rxn_Morgan2FP_separately(
            reactant_smiles, target, rxnfpsize=rxnfpsize, pfpsize=pfpsize, useFeatures=useFeatures
        )
        pfp = np.asarray(pfp, dtype='float32')
        rfp = np.asarray(rfp, dtype='float32')
        rxnfp = pfp - rfp
        return [{
            'input_1': pfp.tolist(),
            'input_2': rxnfp.tolist()
        }]
    
    def transform_output(self, pred):
        """Transforms the output of the prediction into the fast filter score

        Args:
            pred (np.array): numpy array of the output of the prediction

        Returns
            float: the first element of the prediction is returned
        """
        return pred[0]


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
    retroTransformer = RetroTransformer(template_prioritizer=None, fast_filter=None)
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
    global retroTransformer
    return retroTransformer.fast_filter(*args, **kwargs)


@shared_task(bind=True)
def reserve_worker_pool(self):
    """Reserves pool of workers.

    Called by a tb_coordinator to reserve this pool of workers to do a tree
    expansion. This is accomplished by changing what queue(s) this pool
    listens to.

    Returns:
        str: Name of the new queue the worker pool listens to.
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
    """Releases this worker pool so it can listen to the original queues.

    Returns:
        True
    """
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
