from __future__ import absolute_import, unicode_literals, print_function
import os
from celery import Celery
from django.conf import settings

# Set the Django settings module
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'askcos_site.settings')

# Note: cannot use guest for authenticating with broker unless on localhost
# 18.172.0.124 is the IP of RMG.MIT.EDU
SERVERHOST = '18.172.0.124' # rmg.mit.edu
SERVERHOST = '18.63.4.47' # askcos.mit.edu
SERVERHOST = 'localhost'
app = Celery('askcos_site', broker='amqp://worker:askcos@{}:5672'.format(SERVERHOST), 
    backend='django-db',
	include=[
        'askcos_site.askcos_celery.treebuilder.worker',
        'askcos_site.askcos_celery.treebuilder.coordinator',
        'askcos_site.askcos_celery.forwardpredictor.worker',
        'askcos_site.askcos_celery.forwardpredictor.coordinator',
        'askcos_site.askcos_celery.contextrecommender.worker',
        'askcos_site.askcos_celery.chiralretro.coordinator',
        'askcos_site.askcos_celery.chiralretro.worker',
    ]
)

# Using a string here means the worker don't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object('django.conf:settings',)# namespace='CELERY')
app.conf.task_queue_max_priority = 20 # necessary for new tb_worker queues to be priority

if __name__ == '__main__':
    app.start()
