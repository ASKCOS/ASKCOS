"""
The role of a treebuilder coordinator is to take a target compound
and build up the retrosynthetic tree by sending individual chemicals
to workers, which each apply the full set of templates. The coordinator
will keep track of the dictionary and ensure unique IDs in addition
to keeping track of the chemical prices using the Pricer module.

The coordinator, finally, returns a set of buyable trees obtained
from an IDDFS.
"""

from __future__ import absolute_import, unicode_literals, print_function
from django.conf import settings
from celery import shared_task
from celery.signals import celeryd_init
from pymongo import MongoClient
from collections import defaultdict
from celery.result import allow_join_result
# NOTE: allow_join_result is only because the treebuilder worker is separate
from celery.exceptions import Terminated
import time
from askcos_site.askcos_celery.treebuilder.tb_worker import get_top_precursors, reserve_worker_pool, unreserve_worker_pool
from rdkit import RDLogger
import makeit.global_config as gc
from makeit.utilities.buyable.pricer import Pricer
from makeit.retrosynthetic.tree_builder import TreeBuilder
from makeit.synthetic.evaluation.evaluator import Evaluator
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'tb_coordinator'


@celeryd_init.connect
def configure_coordinator(options={}, **kwargs):
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A TREE BUILDER COORDINATOR ###')

    global treeBuilder
    global evaluator

    evaluator = Evaluator(celery=True)

    # Prices
    print('Loading prices...')
    pricer = Pricer()
    pricer.load()
    print('Loaded known prices')
    treeBuilder = TreeBuilder(celery=True, pricer=pricer)

    print('Finished initializing treebuilder coordinator')


@shared_task(bind=True)
def get_buyable_paths(self, smiles, template_prioritization, precursor_prioritization, mincount=0, max_branching=20,
                      max_depth=3, max_ppg=1e8, max_time=60, max_trees=25, reporting_freq=5, known_bad_reactions=[],
                      return_d1_if_no_trees=False, chiral=True, template_count=1e9, forbidden_molecules=[],
                      precursor_score_mode=gc.max, max_cum_template_prob=1, max_natom_dict=defaultdict(lambda: 1e9, {'logic': None}),
                      min_chemical_history_dict={
                          'as_reactant': 1e9, 'as_product': 1e9, 'logic': None},
                      apply_fast_filter=False, filter_threshold=0.8):
    """Get a set of buyable trees for a target compound.

    This function updates its state to provide information about the current
    status of the expansion. Search "on_message task progress" for examples.

    Args:
    smiles (str): SMILES string of the target molecule.
    template_prioritization (str): Keyword for which prioritization method for
        the templates should be used. Keywords can be found in global_config.
    precursor_prioritizer (str): Keyword for which prioritization method for the
        precursors should be used. Keywords can be found in global_config.
    mincount (int, optional): Minimum template popularity. (default: {0})
    max_branching (int, optional): Maximum number of precursor sets to return,
        prioritized using heuristic chemical scoring function. (default: {20})
    max_depth (int, optional): Maximum depth to build the tree out to.
        (default: {3})
    max_ppg (float, optional): Maximum price to consider something buyable.
        (default: {1e8})
    max_time (int, optional): Time (in seconds) for expansion. (default: {60})
    max_trees (int, optional): Maximum number of buyable trees to return. Does
        not affect expansion time. (default: {25})
    reporting_freq (int, optional): Interval (in seconds) for reporting status.
        (default: {5})
    known_bad_reactions (list, optional): Reactant smiles for reactions which we
         already know not to work, so we shouldn't propose them. (default: {[]})
    return_d1_if_no_trees (bool, optional): Unused (default: {False})
    chiral (bool, optional): Whether or not to pay close attention to chirality.
        When False, even achiral templates can lead to accidental inversion of
        chirality in non-reacting parts of the molecule. It is highly
        recommended to keep this as True. (default: {True})
    template_count (int, optional): Maximum number of templates to apply at each
        expansion. (default: {1e9})
    forbidden_molecules (list, optional): Molecules to forbid during expansion,
        represented as a list of SMILES strings. Each string must be
        canonicalized without atom mapping. Forbidden molecules will not be
        allowed as intermediates or leaf nodes. (default: {[]})
    precursor_score_mode (str, optional): Mode to use for precursor scoring when
        using the SCScore prioritizer and multiple reactant fragments must be
        scored together. (default: {gc.max})
    max_cum_template_prob (float, optional): Maximum cumulative template
        probability (i.e., relevance score), used as part of the relevance
        template_prioritizer. (default: {1})
    max_natom_dict (defaultdict, optional): Dictionary defining a potential
        chemical property stopping criterion based on the number of atoms of
        C, N, O, and H in a molecule. The 'logic' keyword of the dict refers to
        how that maximum number of atom information is combined with the
        requirement that chemicals be cheaper than max_ppg.
        (default: {defaultdict(lambda: 1e9, {'logic': None})})
    min_chemical_history_dict (dict, optional): Dictionary defining a potential
        chemical stopping criterion based on the number of times a molecule has
        been seen previously. Always uses logical OR.
        (default: {{'as_reactant':1e9, 'as_product':1e9,'logic':None}})
    apply_fast_filter (bool, optional): Whether to apply the fast filter.
        (default: {False})
    filter_threshold (float, optional): Threshold to use for the fast filter.
        (default: {0.8})

    Returns:
        tree_status ((int, int, dict)): Result of tree_status().
        trees (list of dict): List of dictionaries, where each dictionary
            defines a synthetic route.
    """

    print('Treebuilder coordinator was asked to expand {}'.format(smiles))
    print('Treebuilder coordinator: mincount {}, max_depth {}, max_branching {}, max_ppg {}, max_time {}, max_trees {}'.format(
        mincount, max_depth, max_branching, max_ppg, max_time, max_trees
    ))
    result = treeBuilder.get_buyable_paths(smiles, max_depth=max_depth, max_branching=max_branching, expansion_time=max_time,
                                           template_prioritization=template_prioritization, precursor_prioritization=precursor_prioritization,
                                           mincount=mincount, chiral=chiral, max_trees=max_trees, max_ppg=max_ppg, known_bad_reactions=known_bad_reactions,
                                           template_count=template_count, forbidden_molecules=forbidden_molecules,
                                           precursor_score_mode=precursor_score_mode, max_cum_template_prob=max_cum_template_prob,
                                           max_natom_dict=max_natom_dict, min_chemical_history_dict=min_chemical_history_dict,
                                           apply_fast_filter=apply_fast_filter, filter_threshold=filter_threshold)
    print('Task completed, returning results.')
    return result


@shared_task
def get_template_options():
    """Returns the available options for template prioritizers."""
    return [gc.popularity, gc.relevance, gc.natural]


@shared_task
def get_precursor_options():
    """Returns the available options for precursor prioritizers."""
    return [gc.heuristic, gc.scscore, gc.natural]


@shared_task
def get_precursor_score_options():
    """Returns the available options for precursor scoring modes."""
    retunr[gc.max, gc.mean, gc.geometric, gc.pow8]
