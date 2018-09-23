from __future__ import absolute_import, unicode_literals, print_function
from kombu import Exchange, Queue

TASK_SERIALIZER = 'json'
RESULT_SERIALIZER = 'json'
ACCEPT_CONTENT = ['json']
TIMEZONE = 'US/Eastern'
ENABLE_UTC = True

RESULT_EXPIRES = 1800 # only keep results for 30 minutes max
CELERY_TASK_RESULT_EXPIRES = 1800 # 30 min
CELERY_RESULT_PERSEISTENT = False 

# Custom task queues - necessary to get priority for tree expansion! (RabbitMQ assumed)
TASK_QUEUES = [
    Queue('tb_c_worker', Exchange('tb_c_worker'), routing_key='tb_c_worker', queue_arguments={'x-max-priority':20}),
    Queue('tb_worker', Exchange('tb_worker'), routing_key='tb_worker', queue_arguments={'x-max-priority': 20}),
    
]
CELERY_QUEUES = TASK_QUEUES

# global max priority setting
task_acks_late = True 
CELERY_ACKS_LATE = True
CELERYD_PREFETCH_MULTIPLIER = 1
worker_prefetch_multiplier = 1

# Task routes (to make sure workers are task-specific)
TASK_ROUTES = {
    'askcos_site.askcos_celery.treebuilder.tb_c_worker.get_top_precursors': {'queue': 'tb_c_worker'},
    'askcos_site.askcos_celery.treebuilder.tb_c_worker.fast_filter_check': {'queue': 'tb_c_worker'},
    'askcos_site.askcos_celery.treebuilder.tb_c_worker.apply_one_template_by_idx': {'queue': 'tb_c_worker'},
    'askcos_site.askcos_celery.treebuilder.tb_c_worker.reserve_worker_pool': {'queue': 'tb_c_worker_reservable'},
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
}
CELERY_ROUTES = TASK_ROUTES
