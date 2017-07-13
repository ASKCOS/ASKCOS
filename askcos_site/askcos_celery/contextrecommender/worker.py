'''
The role of a context_worker is to take in an attempted reaction and return
a set of conditions to (try to) run the reaction in. Each worker will
load a pre-trained nearest neighbor model. For each request, this worker
must query the database to get details about the instance.
'''

from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init

NN_PREDICTOR = None 
CORRESPONDING_QUEUE = 'context_worker'

@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    if 'queues' not in options: 
        return 
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A CONTEXT RECOMMENDER WORKER ###')

    global NN_PREDICTOR

    # Get Django settings
    from django.conf import settings

    # Database
    from database import db_client
    db = db_client[settings.INSTANCES['database']]
    INSTANCE_DB = db[settings.INSTANCES['collection']]

    # Setting logging low
    from rdkit import RDLogger
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    # Import nearest-neighbor class
    # from sklearn.externals import joblib
    # from askcos_site.functions.nnPredictor import NNConditionPredictor
    try:
        print('Importing NNConditionPredictor')
        import makeit.webapp.contextmodel as contextmodel
        print('Initializing object')
        NN_PREDICTOR = contextmodel.NNConditionPredictor()
        print('Created NNPredictor object, going to load now...')
        NN_PREDICTOR.load_db_model(db_client[settings.CONTEXT_REC['database']], 
            settings.CONTEXT_REC['model_dir'])
        print('Loaded context recommendation model')
    except Exception as e:
        print(e)
    # # Load all the instance IDs from the test model
    # rxd_ids = []
    # rxn_ids = []
    # with open(settings.CONTEXT_REC['info_path'], 'r') as infile:
    #     rxn_ids.append(infile.readlines()[1:])  # a list of str(rxn_ids) with '\n'
    # for id in rxn_ids[0]:
    #     rxd_ids.append(id.replace('\n', ''))
    # print('Read context recommendation info file from {}'.format(settings.CONTEXT_REC['info_path']))
    # lshf_nn = joblib.load(settings.CONTEXT_REC['model_path'])
    # print('Loaded context recommendation nearest-neighbor model from {}'.format(settings.CONTEXT_REC['model_path']))
    # NN_PREDICTOR = NNConditionPredictor(nn_model=lshf_nn, rxn_ids=rxd_ids, INSTANCE_DB=INSTANCE_DB)
    # NN_PREDICTOR.outputString = False

    print('Context recommendation worker finished loading NN_PREDICTOR')


@shared_task
def get_context_recommendation(rxn, n=1):
    '''Retrieve a context recommendation given the reaction to attempt.

    rxn = [reacants, products], where each is a list of SMILES
    n = number of contexts to return'''

    global NN_PREDICTOR

    print('Context recommender worker got a request for rxn {} and n {}'.format(
        rxn, n))

    return NN_PREDICTOR.step_n_conditions(n, rxn)