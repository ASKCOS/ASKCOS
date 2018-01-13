from __future__ import absolute_import, unicode_literals, print_function
from django.conf import settings
from celery import shared_task
from celery.signals import celeryd_init
from pymongo import MongoClient
from celery.result import allow_join_result
from celery.exceptions import Terminated
import time
from rdkit import RDLogger
from makeit.synthetic.evaluation.tree_evaluator import TreeEvaluator
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'te_coordinator'


@celeryd_init.connect
def configure_coordinator(options={}, **kwargs):
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A TREE EVALUATION COORDINATOR ###')

    global evaluator

    evaluator = TreeEvaluator(celery=True)
    print('### TREE EVALUATION COORDINATOR STARTED UP ###')


@shared_task(bind=True)
def evaluate_tree(self, tree, context_scoring_method='', context_recommender='', forward_scoring_method='', tree_scoring_method='',
                  rank_threshold=5, prob_threshold=0.2, mincount=25, batch_size=500, number_contexts=10, reset=True, template_count = 10000):
    print('Tree evaluation coordinator was asked to evaluate a tree with {} reactions'.format(len(tree)))
    result = evaluator.evaluate_tree(tree,
                                     context_scoring_method=context_scoring_method,
                                     context_recommender=context_recommender,
                                     forward_scoring_method=forward_scoring_method,
                                     tree_scoring_method=tree_scoring_method,
                                     rank_threshold=rank_threshold,
                                     prob_threshold=prob_threshold,
                                     is_target=True,
                                     mincount=mincount,
                                     batch_size=batch_size,
                                     n=number_contexts,
                                     reset=reset,
                                     template_count = template_count)
    print('Task completed, returning results.')
    return result


@shared_task(bind=True)
def evaluate_trees(self, tree_list, context_scoring_method='', context_recommender='', forward_scoring_method='', tree_scoring_method='',
                   rank_threshold=5, prob_threshold=0.2, mincount=25, nproc=1, batch_size=500, n=10, template_count = 10000):
    print('Tree evaluation coordinator was asked to evaluate a list of {} trees.'.format(
        len(tree_list)))
    results = evaluator.evaluate_trees(tree_list,
                                       context_scoring_method=context_scoring_method,
                                       context_recommender=context_recommender,
                                       forward_scoring_method=forward_scoring_method,
                                       tree_scoring_method=tree_scoring_method,
                                       rank_threshold=rank_threshold,
                                       prob_threshold=prob_threshold,
                                       mincount=mincount,
                                       batch_size=batch_size,
                                       n=n,
                                       parallel=False,
                                       template_count = template_count)
    print('Task completed, returning results.')
    return results
