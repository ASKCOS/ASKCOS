from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger
import time
from functools import partial
from rdkit import Chem

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
    # Import as needed
    # from askcos.askcos_site.askcos_celery.impurity.impurity_predictor_worker import predict_reaction
    # from askcos.askcos_site.askcos_celery.impurity.impurity_inspector_worker import inspect_reaction
    # from askcos.askcos_site.askcos_celery.atom_mapper.atom_mapping_worker import get_atom_mapping

    # try:
    #     predictor = predict_reaction
    #     inspector = inspect_reaction
    #     mapper = get_atom_mapping
    # except Exception as e:
    #     print(e)
    #     raise (e)
    print('Initialized')


@shared_task(bind=True)
def get_impurities(self, reactants, reagents='', products='', solvents='',
                   predictor_selection='WLN forward predictor',
                   inspector_selection='Reaxys inspector',
                   mapper_selection='WLN atom mapper',
                   top_k=3, threshold=0.75, check_mapping=True):

    def predictor(reactants_smiles, model=predictor_selection):
        from askcos_site.askcos_celery.impurity.impurity_predictor_worker import predict_reaction
        # print('predictor in impurity_worker', reactants_smiles)
        result = predict_reaction.delay(reactants_smiles, predictor=model)
        return result.get(10)

    def inspector(rxnsmiles, model=inspector_selection):
        from askcos_site.askcos_celery.treebuilder.tb_c_worker import fast_filter_check
        # print('inspector in impurity_worker', rxnsmiles)
        react, prod = rxnsmiles.split('>>')
        result = fast_filter_check.delay(react, prod)
        return result.get(3)

    def mapper(rxnsmiles, model=mapper_selection):
        from askcos_site.askcos_celery.atom_mapper.atom_mapping_worker import get_atom_mapping
        result = get_atom_mapping.delay(rxnsmiles, mapper=model)
        return result.get(10)

    # configure predictor, inspector, and mapper
    # predictor = partial(predict_reaction, predictor=predictor_selection)
    # inspector = partial(inspect_reaction, inspector=inspector_selection)
    # mapper = partial(get_atom_mapping, mapper=mapper_selection)
    # configure impurity predictor

    from makeit.synthetic.impurity.impurity_predictor import ImpurityPredictor
    impurity_predictor = ImpurityPredictor(predictor, inspector, mapper,
                                           topn_outcome=top_k, insp_threshold=threshold,
                                           celery_task=self, check_mapping=check_mapping)
    # make prediction
    outcome = impurity_predictor.predict(reactants, reagents=reagents, products=products, solvents=solvents)
    # return impurity_predictor.predict(reactants, reagents=reagents, products=products, solvents=solvents)
    return outcome

