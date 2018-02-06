import os
time_zero = 0
# Should we use stereochemistry? [should deprecate this setting]

USE_STEREOCHEMISTRY = True

#Output debugging statements
DEBUG = False
#Module options:
#For pathway scoring:
forwardonly = 'Forward only'
templateonly = 'Template only'
product = 'Product'


#For prioritization
all = 'All'
# for precursors
heuristic = 'Heuristic'
scscore = 'SCScore'
mincost = 'MinCost'
mean = 'Mean'
geometric = 'Geometric'
pow8 ='Power of 8'
max = 'Maximum'
# for templates
popularity = 'Popularity'
relevance = 'Relevance'

natural = 'Natural'
# for contexts
probability = 'Probability'
rank = 'Rank'

#For context recommendation
nearest_neighbor = 'Nearest_Neighbor'
neural_network = 'Neural_Network'
#For transformations
template = 'Template'
network = 'Neural_Network'
#Reaction evaluation
fastfilter = 'Fast_Filter'
templatefree = 'Template_Free'
templatebased = 'Template_Based'
forward_scoring_needs_context = {
    'Fast_Filter': False,
    'Template_Free': False,
    'Template_Based': True,
}
forward_scoring_needs_context_necessary_reagent = {
    'Fast_Filter': False,
    'Template_Free': True,
    'Template_Based': True,
}

#Set which modules should be used:
context_module = nearest_neighbor
synth_enumeration = template
retro_enumeration = template
prioritizaton = heuristic
forward_scoring = network

#Use highest protocol in pickle
protocol = -1
data_path = os.path.join(os.path.dirname(__file__),'data')

fingerprint_bits = 256
reaction_fingerprint_bits = 2048

historian_data = os.path.join(data_path, 'historian', 'chemicals.pickle')
reactionhistorian_data = os.path.join(data_path, 'historian', 'reactions.pickle')
pricer_data = os.path.join(data_path,'buyable')
retro_template_data = os.path.join(data_path,'retrosynthetic')
synth_template_data = os.path.join(data_path,'synthetic')
prioritization_data = os.path.join(data_path, 'prioritization')

database = 'reaxys_v2'

#Required database names
MONGO = {
        'path': 'mongodb://guest:guest@askcos2.mit.edu/admin',
        'path_yield_flow':'mongodb://guest:guest@rmg.mit.edu/admin',
        'id': 27017,
        'connect': False
        }


YIELDS = {
    'database': database,
    'collection' : 'yields',
    'data_loc': os.path.join(data_path, 'forward_scoring/yield_est_data'),
    'model_loc': os.path.join(data_path, 'forward_scoring/yield_est_model.h5')
    }

RETRO_TRANSFORMS = {
    'database': database,
    'collection': 'transforms_retro_v8', 
    }
RETRO_TRANSFORMS_CHIRAL = {
    'database': database,
    'collection': 'transforms_retro_v9',
    'mincount': 25,
    'mincount_chiral': 10
}

SYNTH_TRANSFORMS = {
    'database': 'reaxys',
    'collection': 'transforms_forward_v1' ,
    }

INSTANCES = {
    'database': database,
    'collection': 'instances',
    }

REACTIONS = {
    'database': database,
    'collection': 'reactions',
    }

CHEMICALS = {
    'database': database,
    'collection': 'chemicals',
    }

BUYABLES = {
    'database': database,
    'collection': 'buyables',
    }

SOLVENTS = {
    'database': 'reaxys',
    'collection': 'solvents',    
    }

PREDICTOR = {
    'trained_model_path': os.path.join(os.path.dirname(__file__), 'data', 'forward_scoring'),
    'info': '01-23-17, model trained on 80k Reaxys examples, validated on 10k, tested on 10k. Nh1_200, Nh2_200, Nh3_200, l2_0, Nc_5000, enh_weight_0d1, context_weight_50, opt_adadelta, batch_5, moreFeatures'
}

#Hard coded mincounts to maintain compatibility of the relevance method
Relevance_Prioritization = {
    'trained_model_path_True': os.path.join(prioritization_data, 'template_relevance_network_weights.pickle'),
    'min_chiral':10,
    'min':25
    }
SCScore_Prioritiaztion = {
    'trained_model_path_1024bool': os.path.join(prioritization_data, 'scscore', 'model_1024bool.pickle'),
    'trained_model_path_2048bool': os.path.join(prioritization_data, 'scscore', 'model_2048bool.pickle'),
    'trained_model_path_1024uint8': os.path.join(prioritization_data, 'scscore', 'model_1024uint8.pickle')}

MinCost_Prioritiaztion = {
    'trained_model_path': os.path.join(prioritization_data, 'mincost', 'model.hdf5')
    }

CONTEXT_REC = {
    'info_path': os.path.join(data_path,'context', 'RxnID_infoFull.txt'),
    'model_path': os.path.join(data_path,'context', 'fp256noFtr_NN10_BT.pickle'),
    'model_dir': data_path,
    'database': database,
}

NEURALNET_CONTEXT_REC = {
    'info_path': os.path.join(data_path,'context', 'NeuralNet_Cont_Model/'),
    'model_path': os.path.join(data_path,'context', 'NeuralNet_Cont_Model', 'model.json'),
    'weights_path': os.path.join(data_path,'context', 'NeuralNet_Cont_Model', 'weights.h5'),
    'database': database,
}

FLOW_CONDITIONS = {
    'database': database,
    'collection' : 'flow_conditions',
    'data_loc': os.path.join(data_path, 'context/flow_condition_data'), 
    'model_loc': os.path.join(data_path, 'context/flow_model.h5')
    }

FLOW_CONDITIONS_50 = {
    'database':database,
    'collection' : 'flow_conditions_50',
    'data_loc': os.path.join(data_path, 'context/flow_condition_data_50'), 
    'model_loc': os.path.join(data_path, 'context/flow_model_50.h5')
    }
FLOW_CONDITIONS2 = {
    'database': database,
    'collection' : 'flow_conditions',
    'data_loc': os.path.join(data_path, 'context/flow_condition_data2'), 
    'model_loc': os.path.join(data_path, 'context/flow_model2.h5')
    }

FLOW_CONDITIONS2_50 = {
    'database':database,
    'collection' : 'flow_conditions_50',
    'data_loc': os.path.join(data_path, 'context/flow_condition_data2_50'), 
    'model_loc': os.path.join(data_path, 'context/flow_model2_50.h5')
    }

FLOW_CONDITIONS3 = {
    'database': database,
    'collection' : 'flow_conditions',
    'data_loc': os.path.join(data_path, 'context/flow_condition_data3'), 
    'model_loc': os.path.join(data_path, 'context/flow_model3.h5')
    }

FLOW_CONDITIONS3_50 = {
    'database':database,
    'collection' : 'flow_conditions_50',
    'data_loc': os.path.join(data_path, 'context/flow_condition_data3_50'), 
    'model_loc': os.path.join(data_path, 'context/flow_model3_50.h5')
    }

FLOW_CONDITIONS4 = {
    'database': database,
    'collection' : 'flow_conditions',
    'raw_data_loc': os.path.join(data_path, 'context/flow_condition_rawdata4'), 
    'data_loc': os.path.join(data_path, 'context/flow_condition_data4'), 
    'model_loc': os.path.join(data_path, 'context/flow_model4.h5')
    }

FLOW_CONDITIONS4_50 = {
    'database':database,
    'collection' : 'flow_conditions_50',
    'raw_data_loc': os.path.join(data_path, 'context/flow_condition_rawdata4_50'), 
    'data_loc': os.path.join(data_path, 'context/flow_condition_data4_50'),
    'model_loc': os.path.join(data_path, 'context/flow_model4_50.h5')
    }
