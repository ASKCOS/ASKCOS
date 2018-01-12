import os
import makeit.global_config as gc
from pymongo import MongoClient
from logging import MyLogger
from makeit.utilities.buyable.pricer import Pricer
from makeit.synthetic.context.nn_context_recommender import NNContextRecommender
from makeit.synthetic.forward_enumeration.forward_transformer import ForwardTransformer
from makeit.retro_synthetic.retro_transformer import RetroTransformer
from makeit.synthetic.forward_evaluation.template_neuralnet_scorer import TemplateNeuralNetScorer
import sys
model_loader_loc = 'model_loader'


def load_all(retro_mincount=25, retro_mincount_chiral=10):
    MyLogger.print_and_log('Loading models...', model_loader_loc)
    databases = load_Databases()
    pricer = load_Pricer(databases['Chemical_Database'], databases[
                         'Buyable_Database'])

    prioritizer = None
    retroTransformer = None
    synthTransfomrer = None
    scorer = None
    contextRecommender = None

    if(gc.retro_enumeration == gc.template):
        retroTransformer = load_Retro_Transformer(
            databases['Retro_Database'], prioritizer,  mincount=retro_mincount)
    else:
        MyLogger.print_and_log(
            'Invalid retro enumeration method specified. Exiting...', model_loader_loc, level=3)

    models = {
        'retro_transformer': retroTransformer,
        'pricer': pricer,
    }

    MyLogger.print_and_log('All models loaded.', model_loader_loc)

    return models


def load_Retro_Transformer(RETRO_DB=None, mincount=25, mincount_chiral=10, chiral=False):
    '''    
    Load the model and databases required for the retro transformer. Returns the retro transformer, ready to run.
    '''
    MyLogger.print_and_log(
        'Loading retro synthetic template database...', model_loader_loc)
    if not RETRO_DB:
        databases = load_Databases()
        if chiral:
            RETRO_DB = databases['Retro_Database_Chiral']
        else:
            RETRO_DB = databases['Retro_Database']
    retroTransformer = RetroTransformer(
        TEMPLATE_DB=RETRO_DB, mincount=mincount, mincount_chiral=mincount_c)
    retroTransformer.load(chiral=chiral)
    MyLogger.print_and_log(
        'Retro synthetic transformer loaded.', model_loader_loc)
    return retroTransformer


def load_Databases():
    '''
    Load the different databases that will be used: Reactions, Instances, Chemicals, Buyables, Solvents, Retro templates and Synthetic templates
    '''

    db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                            'id'], connect=gc.MONGO['connect'])

    MyLogger.print_and_log('Loading databases...', model_loader_loc)
    db = db_client[gc.REACTIONS['database']]
    REACTION_DB = db[gc.REACTIONS['collection']]

    db = db_client[gc.INSTANCES['database']]
    INSTANCE_DB = db[gc.INSTANCES['collection']]
    db = db_client[gc.CHEMICALS['database']]
    CHEMICAL_DB = db[gc.CHEMICALS['collection']]

    db = db_client[gc.BUYABLES['database']]
    BUYABLE_DB = db[gc.BUYABLES['collection']]
    db = db_client[gc.SOLVENTS['database']]
    SOLVENT_DB = db[gc.SOLVENTS['collection']]

    db = db_client[gc.RETRO_TRANSFORMS['database']]
    RETRO_DB = db[gc.RETRO_TRANSFORMS['collection']]
    RETRO_DB_CHIRAL = db[gc.RETRO_TRANSFORMS_CHIRAL['collection']]
    db = db_client[gc.SYNTH_TRANSFORMS['database']]
    SYNTH_DB = db[gc.SYNTH_TRANSFORMS['collection']]

    databases = {
        'Reaction_Database': REACTION_DB,
        'Instance_Database': INSTANCE_DB,
        'Chemical_Database': CHEMICAL_DB,
        'Buyable_Database': BUYABLE_DB,
        'Solvent_Database': SOLVENT_DB,
        'Retro_Database': RETRO_DB,
        'Retro_Database_Chiral': RETRO_DB_CHIRAL,
        'Synth_Database': SYNTH_DB
    }
    MyLogger.print_and_log('Databases loaded.', model_loader_loc)
    return databases


def load_Pricer(chemical_database, buyable_database):
    '''
    Load a pricer using the chemicals database and database of buyable chemicals
    '''
    MyLogger.print_and_log('Loading pricing model...', model_loader_loc)
    pricerModel = Pricer()
    pricerModel.load(chemical_database, buyable_database)
    MyLogger.print_and_log('Pricer Loaded.', model_loader_loc)
    return pricerModel


def load_Forward_Transformer(SYNTH_DB, mincount=100):
    '''
    Load the forward prediction neural network
    '''
    MyLogger.print_and_log(
        'Loading forward prediction model...', model_loader_loc)
    transformer = ForwardTransformer(TEMPLATE_DB=SYNTH_DB, mincount=mincount)
    transformer.load()
    MyLogger.print_and_log('Forward transformer loaded.', model_loader_loc)
    return transformer


def load_fastfilter():
    # Still has to be implemented
    return None


def load_templatebased(mincount=25, celery=False):
    transformer = None
    databases = load_Databases()
    if not celery:
        transformer = load_Forward_Transformer(
            databases['Synth_Database'], mincount=mincount)

    scorer = TemplateNeuralNetScorer(
        forward_transformer=transformer, celery=celery)
    scorer.load(databases['Solvent_Database'],
                gc.PREDICTOR['trained_model_path'])
    return scorer


def load_templatefree():
    # Still has to be implemented
    return None


def load_Context_Recommender(context_recommender, max_contexts=10):
    '''
    Load the context recommendation model
    '''
    MyLogger.print_and_log('Loading context recommendation model: {}...'.format(
        context_recommender), model_loader_loc)
    if context_recommender == gc.nearest_neighbor:
        recommender = NNContextRecommender(max_contexts=max_contexts)
        recommender.load(model_path=gc.CONTEXT_REC[
                         'model_path'], info_path=gc.CONTEXT_REC['info_path'])
    else:
        raise NotImplementedError
    MyLogger.print_and_log('Context recommender loaded.', model_loader_loc)
    return recommender
