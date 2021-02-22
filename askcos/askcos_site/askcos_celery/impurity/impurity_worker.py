from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger

from makeit.synthetic.impurity.impurity_predictor import ImpurityPredictor
from ..atom_mapper.atom_mapping_worker import get_atom_mapping
from ..impurity.impurity_predictor_worker import predict_reaction
from ..treebuilder.tb_c_worker import fast_filter_check

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'impurity_worker'


@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    print(options)
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A IMPURITY PREDICTOR WORKER ###')
    print('Initialized')


@shared_task(bind=True)
def get_impurities(self, reactants, reagents='', products='', solvents='',
                   predictor_selection='WLN forward predictor',
                   inspector_selection='Reaxys inspector',
                   mapper_selection='WLN atom mapper',
                   top_k=3, threshold=0.75, check_mapping=True):

    def predictor(reactants_smiles, model=predictor_selection):
        result = predict_reaction.delay(reactants_smiles, predictor=model)
        return result.get(10)

    def inspector(rxnsmiles, model=inspector_selection):
        if model == 'Reaxys inspector':
            react, prod = rxnsmiles.split('>>')
            result = fast_filter_check.delay(react, prod)
        else:
            raise NotImplementedError('{0} is not yet supported for impurity prediction.'.format(model))
        return result.get(3)

    def mapper(rxnsmiles, model=mapper_selection):
        result = get_atom_mapping.delay(rxnsmiles, mapper=model)
        return result.get(10)

    impurity_predictor = ImpurityPredictor(predictor, inspector, mapper,
                                           topn_outcome=top_k, insp_threshold=threshold,
                                           celery_task=self, check_mapping=check_mapping)
    # make prediction
    return impurity_predictor.predict(reactants, reagents=reagents, products=products, solvents=solvents)
