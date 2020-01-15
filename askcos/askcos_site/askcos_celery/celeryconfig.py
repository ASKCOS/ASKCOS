from __future__ import absolute_import, unicode_literals, print_function
from kombu import Exchange, Queue

CELERY_TASK_SERIALIZER = 'json'
CELERY_RESULT_SERIALIZER = 'json'
CELERY_ACCEPT_CONTENT = ['json']
CELERY_TIMEZONE = 'US/Eastern'
CELERY_ENABLE_UTC = True

CELERY_TASK_RESULT_EXPIRES = 1800 # 30 min
CELERY_RESULT_PERSISTENT = False

CELERY_TASK_QUEUES = [
    Queue('tb_c_worker', Exchange('tb_c_worker'), routing_key='tb_c_worker', queue_arguments={'x-max-priority':20}),
    Queue('tb_worker', Exchange('tb_worker'), routing_key='tb_worker', queue_arguments={'x-max-priority': 20}),

]


CELERY_BROKER_HEARTBEAT = 0

# global max priority setting
CELERY_TASK_ACKS_LATE = True
CELERYD_PREFETCH_MULTIPLIER = 1

# Task routes (to make sure workers are task-specific)
CELERY_TASK_ROUTES = {
    'askcos_site.askcos_celery.treebuilder.tb_c_worker_preload.get_top_precursors': {'queue': 'tb_c_worker_preload'},
    'askcos_site.askcos_celery.treebuilder.tb_c_worker.get_top_precursors': {'queue': 'tb_c_worker'},
    'askcos_site.askcos_celery.treebuilder.tb_c_worker.fast_filter_check': {'queue': 'tb_c_worker'},
    'askcos_site.askcos_celery.treebuilder.tb_c_worker.apply_one_template_by_idx': {'queue': 'tb_c_worker'},
    'askcos_site.askcos_celery.treebuilder.tb_c_worker.reserve_worker_pool': {'queue': 'tb_c_worker_reservable'},
    'askcos_site.askcos_celery.treebuilder.tb_c_worker.template_relevance': {'queue': 'tb_c_worker'},
    'askcos_site.askcos_celery.treebuilder.tb_worker.get_top_precursors': {'queue': 'tb_worker'},
    'askcos_site.askcos_celery.treebuilder.tb_worker.reserve_worker_pool': {'queue': 'tb_worker_reservable'},
    'askcos_site.askcos_celery.treebuilder.tb_coordinator.*':{'queue': 'tb_coordinator'},
    'askcos_site.askcos_celery.treebuilder.tb_coordinator_mcts.*':{'queue': 'tb_coordinator_mcts'},
    'askcos_site.askcos_celery.treeevaluator.tree_evaluation_coordinator.*':{'queue':'te_coordinator'},
    'askcos_site.askcos_celery.treeevaluator.scoring_coordinator.*':{'queue':'sc_coordinator'},
    'askcos_site.askcos_celery.treeevaluator.forward_trans_worker.*': {'queue':'ft_worker'},
    'askcos_site.askcos_celery.contextrecommender.cr_nn_worker.*': {'queue':'cr_nn_worker'},
    'askcos_site.askcos_celery.contextrecommender.cr_network_worker.*': {'queue':'cr_network_worker'},
    'askcos_site.askcos_celery.contextrecommender.cr_coordinator.*':{'queue':'cr_coordinator'},
    'askcos_site.askcos_celery.siteselectivity.sites_worker.*':{'queue':'sites_worker'},
    'askcos_site.askcos_celery.impurity.impurity_worker.*':{'queue':'impurity_worker'},
    'askcos_site.askcos_celery.atom_mapper.atom_mapping_worker.*':{'queue':'atom_mapping_worker'},
    'askcos_site.askcos_celery.impurity.impurity_inspector_worker.*': {'queue': 'atom_mapping_worker'},
    'askcos_site.askcos_celery.impurity.impurity_predictor_worker.*': {'queue': 'atom_mapping_worker'},

}
