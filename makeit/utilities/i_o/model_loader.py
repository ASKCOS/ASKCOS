import os
import global_config as gc
from pymongo import MongoClient
from logging import MyLogger
import models.transformer as transformer
from models.pricer import Pricer
from models.forwardPredictor import ForwardPredictor
from models.contextmodel import NNConditionPredictor as ConditionPredictor
from models.scorer import Scorer
from multiprocessing import Process
#from application.MyQueue import Queue
from multiprocessing import JoinableQueue as Queue
from functools import partial
import threading
import sys
model_loader_loc = 'model_loader'

"""
def load_all(retro_mincount,synth_mincount):
    #parallel does not work with MongoClient objects
    output = Queue()
    
    #load databases first: fast, and databases are required to load the models.
    databases = load_Databases()
    
    #set up parallel processes for loading models
    retro_process = Process(target = load_Retro_Transformer, args = (databases['Retro_Database'], output), kwargs = {'mincount_retro' : retro_mincount})
    synth_process = Process(target = load_Forward_Transformer, args = (databases['Synth_Database'], output), kwargs = {'mincount_synth' : synth_mincount})
    pricer_process = Process(target = load_Pricer, args = (databases['Chemical_Database'], databases['Buyable_Database'], output))
    #pred_process = Process(target = load_Forward_Predictor, args = (databases['Synth_Database'], databases['Solvent_Database'], output), kwargs = {'mincount_synth' : synth_mincount})
    #processes = {retro_process, synth_process, pricer_process, pred_process}
    processes = {retro_process, synth_process,pricer_process}
    #Run all processes
    for p in processes:
        p.start()
    # Exit the completed processes
     
    for p in processes:
        p.join()

    models = [output.get() for p in processes]
    
    print(models)

"""
   
def load_Retro_Transformer(RETRO_DB, mincount_retro = 250, nb_workers = 1):
    '''    
    Load the model and databases required for the retro transformer. Returns the retro transformer, ready to run.
    '''
    MyLogger.print_and_log('Loading retrosynthetic template database...',model_loader_loc)
    retroTransformer = transformer.Transformer(
        parallel = False if nb_workers == 1 else True, 
        nb_workers = nb_workers,
        )
    retroTransformer.load(RETRO_DB, mincount = mincount_retro, get_retro = True, get_synth = False)
    MyLogger.print_and_log('Retrosynthetic transformer loaded.',model_loader_loc)
    #output.put(retroTransformer)
    return retroTransformer
    
def load_Forward_Transformer(SYNTH_DB, mincount_synth = 100):
    '''
    Load the model and databases required for the forward transformer. Returns forward transformer, ready to run.
    '''
    MyLogger.print_and_log('Loading synthetic template database...',model_loader_loc)
    synth_Transformer = transformer.Transformer()
    synth_Transformer.load(SYNTH_DB, mincount_synth, get_retro = False, get_synth = True)
    MyLogger.print_and_log('Synthetic transformer loaded.', model_loader_loc)
    #output.put(synth_Transformer)
    return synth_Transformer
        
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
    #output.put(pricerModel) 
    return pricerModel

def load_Forward_Predictor(SYNTH_DB, mincount_synth = 100):
    '''
    Load the forward prediction neural network
    '''
    MyLogger.print_and_log('Loading forward prediction model...',model_loader_loc)
    predictor = ForwardPredictor(TRANSFORM_DB = SYNTH_DB, mincount = mincount_synth)
    predictor.load_templates()
    MyLogger.print_and_log('Forward predictor loaded.',model_loader_loc)
    #output.put(predictor)
    return predictor

def load_Scorer():
    '''
    Load the neural network for scoring a reaction
    '''
    MyLogger.print_and_log('Loading reaction scoring neural network...',model_loader_loc)
    scorer = Scorer()
    scorer.load_model(gc.PREDICTOR['trained_model_path'])
    
    MyLogger.print_and_log('Reaction scoring neural network loaded.',model_loader_loc)
    
    return scorer
    
def load_Context_Recommender(REACTION_DB, INSTANCE_DB, CHEMICAL_DB, SOLVENT_DB, max_total_contexts):
    '''
    Load the context recommendation model
    '''
    
    MyLogger.print_and_log('Loading context recommendation model...', model_loader_loc)
    recommender = ConditionPredictor(max_total_contexts = max_total_contexts)
    recommender.load_db_model(gc.CONTEXT_REC['model_path'], gc.CONTEXT_REC['info_path'], REACTION_DB, INSTANCE_DB, CHEMICAL_DB, SOLVENT_DB)
    MyLogger.print_and_log('Context recommender loaded.', model_loader_loc)
    #output.put(recommender)
    return recommender