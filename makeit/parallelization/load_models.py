import global_config as gc
import time
from multiprocessing import Process, Manager
from synthetic.forward_enumeration.forward_transformer import ForwardTransformer
from retro_synthetic.retro_transformer import RetroTransformer
from utilities.buyable.pricer import Pricer
from synthetic.forward_evaluation.scorer import Scorer
from synthetic.context.nn_context_recommender import NNConditionPredictor
from utilities.i_o.logging import MyLogger
load_models_loc = 'load_models'

def load_all(nproc = 3, mincount_f = 4, mincount_r = 4, NN_model_path = "", context_model_path = "", context_info_path = "", max_ppg = 1e10, max_total_contexts = 10):
    '''
    Load all models in one go and in parallel. Optimal number of processors: 3.
    '''
    
    queue = Manager().Queue()
    
    forward_done = Manager().Value('i', 0)
    forward_transformer = ForwardTransformer(mincount = mincount_f, done = forward_done)
    f_trans = Process(target = forward_transformer.load, kwargs = {'queue':queue})
    
    retro_done = Manager().Value('i', 0)
    retro_transformer = RetroTransformer(mincount = mincount_r, done = retro_done)
    r_trans = Process(target = retro_transformer.load, kwargs = {'queue':queue})
    
    scorer_done = Manager().Value('i', 0)
    scorer = Scorer(done = scorer_done)
    scor = Process(target = scorer.load, kwargs = {'folder': NN_model_path, 'queue':queue})
    
    pricer_done = Manager().Value('i', 0)
    pricer = Pricer(max_ppg = max_ppg, done = pricer_done)
    pric = Process(target = pricer.load, kwargs = {'queue':queue})
    
    context_done = Manager().Value('i', 0)
    context_recommender = NNConditionPredictor(max_total_contexts = max_total_contexts, done = context_done)
    cont_rec = Process(target = context_recommender.load, kwargs = {'model_path': context_model_path, 'info_path': context_info_path, 'queue':queue})
    
    processes = [r_trans, cont_rec, f_trans, pric, scor]
    
    checks = [retro_done.value, context_done.value, forward_done.value, pricer_done.value, scorer_done.value]
    last_start = 0
    #for each available node: start up one of the processes. If the number of nodes specified > number of models to be loaded,
    #five nodes will be used. Otherwise, the max number of nodes will be started up and proceed to next block.
    for i,process in enumerate(processes):
        if(i < nproc):
            process.start()
            last_start = i
        else:
            break
    
    #As long as not all modules loaded: look for modules that have finished and thus freed up a node. If a process finishes,
    #Start up the next process in line, until all processes have been started. Then just check whether all modules are loaded.
    while not all(checks):
        checks = [retro_done.value, context_done.value, forward_done.value, pricer_done.value, scorer_done.value]
        time.sleep(1)
        for i,check in enumerate(checks):
            if(i < last_start):
                if check and (last_start+1 < len(processes)):
                    #end finished process.
                    processes[i].join()
                    processes[i].terminate()
                    #start new process.
                    processes[last_start+1].start()
                    last_start += 1
    
    models = {}
    while not queue.empty():
        (tag, model) = queue.get(0.1)
        models[tag] = model
    
    models = {
        'retro_transformer':retro_transformer,
        'forward_transformer':forward_transformer,
        'scorer':scorer,
        'pricer':pricer,
        'context_recommender':context_recommender
        }
      
    #clean up all running processes.
    for process in processes:
        process.join()
        process.terminate()
 
    return models

if __name__ == '__main__':
    MyLogger.initialize_logFile()
    load_all(nproc = 4, mincount_f = 10, mincount_r = 10, NN_model_path = gc.PREDICTOR['trained_model_path'], 
            context_model_path = gc.CONTEXT_REC['model_path'], context_info_path = gc.CONTEXT_REC['info_path'])
