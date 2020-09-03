from __future__ import absolute_import, unicode_literals, print_function
import os
from celery import Celery
from django.conf import settings

# Set the Django settings module
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'askcos_site.settings')

# Note: cannot use guest for authenticating with broker unless on localhost
REDIS_HOST = os.environ.get('REDIS_HOST', 'localhost')
RABBIT_HOST = os.environ.get('RABBIT_HOST', 'localhost')
app = Celery('askcos_site', broker='amqp://{}:5672'.format(RABBIT_HOST),
    backend='redis://{}:6379'.format(REDIS_HOST),
    include=[
        'askcos_site.askcos_celery.treebuilder.tb_worker',
        'askcos_site.askcos_celery.treebuilder.tb_c_worker',
        'askcos_site.askcos_celery.treebuilder.tb_coordinator',
        'askcos_site.askcos_celery.treebuilder.tb_coordinator_mcts',
        'askcos_site.askcos_celery.contextrecommender.cr_coordinator',
        'askcos_site.askcos_celery.contextrecommender.cr_nn_worker',
        'askcos_site.askcos_celery.contextrecommender.cr_network_worker',
        'askcos_site.askcos_celery.treeevaluator.forward_trans_worker',
        'askcos_site.askcos_celery.treeevaluator.scoring_coordinator',
        'askcos_site.askcos_celery.treeevaluator.tree_evaluation_coordinator',
    ]
)

# Using a string here means the worker don't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object('django.conf:settings', namespace='CELERY')
app.conf.task_queue_max_priority = 20 # necessary for new tb_worker queues to be priority

if __name__ == '__main__':
    app.start()
