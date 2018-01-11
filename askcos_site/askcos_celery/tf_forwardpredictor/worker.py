'''
Template-free forward predictor
'''

from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init

tffp = None 
CORRESPONDING_QUEUE = 'tffp_worker'

@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    if 'queues' not in options: 
        return 
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A TEMPLATE-FREE FORWARD PREDICTOR WORKER ###')

    global tffp

    # Setting logging low
    from rdkit import RDLogger
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    # Import as needed
    from makeit.synthetic.forward_evaluation.rexgen_release.predict import TFFP 
    print('Imported TFFP')
    try:
        tffp = TFFP()
    except Exception as e:
        print(e)
        raise(e)
    print('Initialized')
    tffp.predict('CCCO.CCCBr')

    print('Finished configuring TFFP worker')

@shared_task
def get_outcomes(reactants, mincount=0, top_n=10):
    global tffp 
    return tffp.predict(reactants, top_n=top_n)
