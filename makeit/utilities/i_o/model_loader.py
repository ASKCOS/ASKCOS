import os
import global_config as gc
from pymongo import MongoClient
from logging import MyLogger
from utilities.buyable.pricer import Pricer
from synthetic.context.nn_context_recommender import NNContextRecommender
from synthetic.forward_evaluation.scorer import Scorer
from synthetic.forward_enumeration.forward_transformer import ForwardTransformer
from retro_synthetic.retro_transformer import RetroTransformer
from retro_synthetic.heuristic_prioritizer import HeuristicPrioritizer


from synthetic.forward_evaluation.template_neuralnet_scorer import TemplateNeuralNetScorer
from multiprocessing import Process
from multiprocessing import JoinableQueue as Queue
from functools import partial
import threading
import sys
model_loader_loc = 'model_loader'


def load_all(retro_mincount = 100, synth_mincount = 100, max_contexts = 10):
        MyLogger.print_and_log('Loading models...', model_loader_loc)
        databases = load_Databases()
        pricer = load_Pricer(databases['Chemical_Database'], databases['Buyable_Database'])
        
        prioritizer = None
        retroTransformer = None
        synthTransfomrer = None
        scorer = None
        contextRecommender = None
        
        if(gc.prioritizaton == gc.heuristic):
            prioritizer = HeuristicPrioritizer()
        else:
            MyLogger.print_and_log('Invalid prioritization method specified. Exiting...', model_loader_loc, level = 3)
        
        if(gc.retro_enumeration == gc.template):
            retroTransformer = load_Retro_Transformer(databases['Retro_Database'],prioritizer,  mincount_retro = retro_mincount)
        else:
            MyLogger.print_and_log('Invalid retro enumeration method specified. Exiting...', model_loader_loc, level = 3)
        
        if(gc.forward_enumeration == gc.template):
            synthTransformer = load_Forward_Transformer(databases['Synth_Database'], mincount_synth = synth_mincount)
        else:
            MyLogger.print_and_log('Invalid forward enumeration method specified. Exiting...', model_loader_loc, level = 3)
        
        if(gc.forward_scoring == gc.network):
            scorer = load_Scorer()
        else:
            MyLogger.print_and_log('Invalid scoring method specified. Exiting...', model_loader_loc, level = 3)
            
        if(gc.context_module == gc.nearest_neighbor):
            contextRecommender = load_Context_Recommender(max_total_contexts = max_contexts)
        else:
            MyLogger.print_and_log('Invalid context recommendation method specified. Exiting...', model_loader_loc, level = 3)
        
        models ={
            'retro_transformer':retroTransformer,
            'prioritizer': prioritizer,
            'synthetic_transformer':synthTransformer,
            'pricer':pricer,
            'context_recommender':contextRecommender,
            'scorer':scorer
            }
        
        MyLogger.print_and_log('All models loaded.', model_loader_loc)
        
        return models

   
def load_Retro_Transformer(RETRO_DB, prioritizer, mincount_retro = 250):
    '''    
    Load the model and databases required for the retro transformer. Returns the retro transformer, ready to run.
    '''
    MyLogger.print_and_log('Loading retrosynthetic template database...',model_loader_loc)
    retroTransformer = RetroTransformer(prioritizer = prioritizer, TEMPLATE_DB = RETRO_DB, mincount = mincount_retro)
    retroTransformer.load()
    MyLogger.print_and_log('Retrosynthetic transformer loaded.',model_loader_loc)
    return retroTransformer
    
def load_Databases():
    '''
    Load the different databases that will be used: Reactions, Instances, Chemicals, Buyables, Solvents, Retro templates and Synthetic templates
    '''
    
    db_client = MongoClient(gc.MONGO['path'], gc.MONGO['id'], connect = gc.MONGO['connect'])
    
    MyLogger.print_and_log('Loading databases...',model_loader_loc)
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
    
    db = db_client[gc.SYNTH_TRANSFORMS['database']]
    SYNTH_DB = db[gc.SYNTH_TRANSFORMS['collection']]
    
    databases = {
        'Reaction_Database' : REACTION_DB,
        'Instance_Database' : INSTANCE_DB,
        'Chemical_Database' : CHEMICAL_DB,
        'Buyable_Database' : BUYABLE_DB,
        'Solvent_Database' : SOLVENT_DB,
        'Retro_Database' : RETRO_DB,
        'Synth_Database' : SYNTH_DB
        }        
    MyLogger.print_and_log('Databases loaded.', model_loader_loc)        
    return databases
    
def load_Pricer(chemical_database, buyable_database):
    '''
    Load a pricer using the chemicals database and database of buyable chemicals
    '''
    MyLogger.print_and_log('Loading pricing model...',model_loader_loc)
    pricerModel = Pricer()
    pricerModel.load(chemical_database, buyable_database)
    MyLogger.print_and_log('Pricer Loaded.',model_loader_loc)
    return pricerModel

def load_Forward_Transformer(SYNTH_DB, mincount_synth = 100):
    '''
    Load the forward prediction neural network
    '''
    MyLogger.print_and_log('Loading forward prediction model...',model_loader_loc)
    transformer = ForwardTransformer(TEMPLATE_DB = SYNTH_DB, mincount = mincount_synth)
    transformer.load()
    MyLogger.print_and_log('Forward transformer loaded.',model_loader_loc)
    return transformer

def load_Scorer():
    '''
    Load the neural network for scoring a reaction
    '''
    MyLogger.print_and_log('Loading reaction scoring neural network...',model_loader_loc)
    scorer = Scorer()
    scorer.load(folder = gc.PREDICTOR['trained_model_path'])
    
    MyLogger.print_and_log('Reaction scoring neural network loaded.',model_loader_loc)
    return scorer

def load_fastfilter():
    #Still has to be implemented
    return None

def load_templatebased(chiral = False, mincount = 25, celery = False):
    transformer = None
    databases = load_Databases()
    if not celery:
        transformer = ForwardTransformer(mincount = mincount)
        transformer.load(chiral = chiral)
        
    scorer = TemplateNeuralNetScorer(forward_transformer = transformer, celery = celery)
    scorer.load(databases['Solvent_Database'], gc.PREDICTOR['trained_model_path'])
    return scorer

def load_templatefree():
    #Still has to be implemented
    return None

def load_Context_Recommender(max_total_contexts):
    '''
    Load the context recommendation model
    '''
    
    MyLogger.print_and_log('Loading context recommendation model...', model_loader_loc)
    recommender = NNContextRecommender(max_total_contexts = max_total_contexts)
    recommender.load(model_path = gc.CONTEXT_REC['model_path'], info_path = gc.CONTEXT_REC['info_path'])
    MyLogger.print_and_log('Context recommender loaded.', model_loader_loc)
    return recommender