from __future__ import absolute_import, unicode_literals, print_function
import os
from celery import Celery
from django.conf import settings

# Set the Django settings module
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'askcos_site.settings')

app = Celery('askcos_site', broker='amqp://localhost', backend='django-db',
	include=[
    'askcos_site.askcos_celery.treebuilder.worker',
    'askcos_site.askcos_celery.treebuilder.coordinator',
    'askcos_site.askcos_celery.forwardpredictor.worker',
    'askcos_site.askcos_celery.forwardpredictor.coordinator',
    'askcos_site.askcos_celery.contextrecommender.worker',
])

# Using a string here means the worker don't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object('django.conf:settings',)# namespace='CELERY')


if __name__ == '__main__':
    app.start()
