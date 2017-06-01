from __future__ import absolute_import, unicode_literals, print_function
from kombu import Exchange, Queue

TASK_SERIALIZER = 'json'
RESULT_SERIALIZER = 'json'
ACCEPT_CONTENT = ['json']
TIMEZONE = 'US/Eastern'
ENABLE_UTC = True

RESULT_EXPIRES = 300 # only keep results for 5 minutes max

# Custom task queues - necessary to get priority for tree expansion! (RabbitMQ assumed)
TASK_QUEUES = [
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
    'askcos_site.askcos_celery.treebuilder.worker.get_top_precursors': {'queue': 'tb_worker'},
    'askcos_site.askcos_celery.treebuilder.worker.reserve_worker_pool': {'queue': 'tb_worker_reservable'},
    'askcos_site.askcos_celery.treebuilder.coordinator.*': {'queue': 'tb_coordinator'},
    'askcos_site.askcos_celery.forwardpredictor.worker.*': {'queue': 'fp_worker'},
    'askcos_site.askcos_celery.forwardpredictor.coordinator.*': {'queue': 'fp_coordinator'}, 
    'askcos_site.askcos_celery.contextrecommender.worker.*': {'queue': 'context_worker'},   
    'askcos_site.askcos_celery.chiralretro.coordinator.*': {'queue': 'cr_coordinator'},
    'askcos_site.askcos_celery.chiralretro.worker.*': {'queue': 'cr_worker'},
}
CELERY_ROUTES = TASK_ROUTES

