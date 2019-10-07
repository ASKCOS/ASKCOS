'''
Site selectivity predictor
'''
from __future__ import absolute_import, unicode_literals, print_function
from django.conf import settings
from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

sites_pred = None
CORRESPONDING_QUEUE = 'sites_worker'

@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    print(options)
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A SITE PREDICTOR WORKER ###')
    global sites_pred
    # Import as needed
    from makeit.synthetic.selectivity.site_selectivity import Site_Predictor
    try:
        sites_pred = Site_Predictor()
    except Exception as e:
        raise(e)
    print('Initialized')
    print(sites_pred.predict('Cc1ccccc1'))
    print('Finished configuring sites worker')


@shared_task
def get_sites(smi):
    global sites_pred
    print('site selectivity got a request {}'.format(smi))
    res = sites_pred.predict(smi)
    return res