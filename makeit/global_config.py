import os
time_zero = 0

### TODO: deprecate this - always use stereochemsitry!
# Should we use stereochemistry?
USE_STEREOCHEMISTRY = True

# Output debugging statements
DEBUG = False

# Whether to preload all templates for the retrotransformer
PRELOAD_TEMPLATES = False

################################################################################
# Options for different modules, defined as strings
################################################################################

# For pathway scoring:
forwardonly = 'Forward only'
templateonly = 'Template only'
product = 'Product'

# For prioritization (precursors and templates)
all = 'All'

# For precursor prioritization
relevanceheuristic = 'RelevanceHeuristic'
heuristic = 'Heuristic'
scscore = 'SCScore'
mincost = 'MinCost'
mean = 'Mean'
geometric = 'Geometric'
pow8 ='Power of 8'
max = 'Maximum'

# For template prioritization
popularity = 'Popularity'
relevance = 'Relevance'
natural = 'Natural'

# For deciding the best context
probability = 'Probability'
rank = 'Rank'

# For context recommendation
nearest_neighbor = 'Nearest_Neighbor'
neural_network = 'Neural_Network'

# For forward prediction
template = 'Template'
network = 'Neural_Network'

# For reaction evaluation
fastfilter = 'Fast_Filter'
templatefree = 'Template_Free'
templatebased = 'Template_Based'
forward_scoring_needs_context = {
    'Fast_Filter': False,
    'Template_Free': True,
    'Template_Based': True,
}
forward_scoring_needs_context_necessary_reagent = {
    'Fast_Filter': False,
    'Template_Free': True,
    'Template_Based': True,
}

# Set which modules should be used as defaults
context_module = nearest_neighbor
synth_enumeration = template
retro_enumeration = template
prioritizaton = heuristic
forward_scoring = network


################################################################################
# Define data file locations
################################################################################

data_path = os.path.join(os.path.dirname(__file__),'data')
local_db_dumps = os.path.join(data_path, 'local_db_dumps')

fingerprint_bits = 256
reaction_fingerprint_bits = 2048

historian_data = os.path.join(data_path, 'historian', 'chemicals.pickle')
reactionhistorian_data = os.path.join(data_path, 'historian', 'reactions.pickle')
retro_template_data = os.path.join(data_path,'retrosynthetic')
synth_template_data = os.path.join(data_path,'synthetic')
prioritization_data = os.path.join(data_path, 'prioritization')

database = 'reaxys_v2'

################################################################################
# Define databases (should be nonessential if all local files present)
################################################################################

MONGO_HOST = os.environ.get('MONGO_HOST', 'MONGO_HOST')
MONGO_USER = os.environ.get('MONGO_USER', 'USERNAME')
MONGO_PW = os.environ.get('MONGO_PW', 'PASSWORD')

# TODO: change this to your local Mongo DB!
MONGO = {
    'path': 'mongodb://{}:{}@{}'.format(MONGO_USER, MONGO_PW, MONGO_HOST),
    'id': 27017,
    'connect': False
}

# TODO: deprecate achiral transforms
RETRO_TRANSFORMS = {
    'database': database,
    'collection': 'transforms_retro_v6',
    'mincount': 25,
}

RETRO_TRANSFORMS_CHIRAL = {
    'file_name': 'retrotransformer_chiral_using_reaxys_v2-transforms_retro_v9_mincount10_mincountchiral5.pkl',
    'database': 'askcos',
    'collection': 'retro_templates',
    'mincount': 10,
    'mincount_chiral': 5
}

SYNTH_TRANSFORMS = {
    'database': 'reaxys',
    'collection': 'transforms_forward_v1' ,
    'mincount': 25,
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

CHEMICAL_HISTORY = {
    'database': database,
    'collection': 'chemical_history'
}

BUYABLES = {
    'file_name': 'pricer_using_reaxys_v2-chemicals_and_reaxys_v2-buyables.pkl',
    'database': 'askcos',
    'collection': 'buyables',
}

SOLVENTS = {
    'database': 'reaxys',
    'collection': 'solvents',
}

# Template-based forward predictor
PREDICTOR = {
    'trained_model_path': os.path.join(os.path.dirname(__file__), 'data', 'forward_scoring'),
    'info': '01-23-17, model trained on 80k Reaxys examples, validated on 10k, tested on 10k. Nh1_200, Nh2_200, Nh3_200, l2_0, Nc_5000, enh_weight_0d1, context_weight_50, opt_adadelta, batch_5, moreFeatures'
}

# Fast filter evaluation
FAST_FILTER_MODEL = {
    'trained_model_path': os.path.join(os.path.dirname(__file__), 'data', 'fast_filter','fast_filter_cleandata.h5'),
}

# Hard coded mincounts to maintain compatibility of the relevance method (weights are numpy matrices)
Relevance_Prioritization = {
    'trained_model_path_True': os.path.join(prioritization_data, 'template_relevance_network_weights_v9_10_5.pickle'),
    'output_size': 163723,
    'min_chiral':5,
    'min':10
}

# Relevance_Prioritization_OLD = {
#     'trained_model_path_True': os.path.join(prioritization_data, 'template_relevance_network_weights_v9_25_10.pickle'),
#     'output_size': 61142,
#     'min_chiral':10,
#     'min':25
# }

# Different SCScore models that are all functionally similary
SCScore_Prioritiaztion = {
    'trained_model_path_1024bool': os.path.join(prioritization_data, 'scscore', 'model_1024bool.pickle'),
    'trained_model_path_2048bool': os.path.join(prioritization_data, 'scscore', 'model_2048bool.pickle'),
    'trained_model_path_1024uint8': os.path.join(prioritization_data, 'scscore', 'model_1024uint8.pickle')}

MinCost_Prioritiaztion = {
    'trained_model_path': os.path.join(prioritization_data, 'mincost', 'model.hdf5')
}

CONTEXT_REC = {
    'info_path': os.path.join(data_path, 'context', 'RxnID_infoFull.txt'),
    'model_path': os.path.join(data_path, 'context', 'fp256noFtr_NN10_BT.pickle'),
    'model_dir': data_path,
    'database': database,
}

NEURALNET_CONTEXT_REC = {
    'info_path': os.path.join(data_path,'context', 'NeuralNet_Cont_Model/'),
    'model_path': os.path.join(data_path,'context', 'NeuralNet_Cont_Model', 'model.json'),
    'weights_path': os.path.join(data_path,'context', 'NeuralNet_Cont_Model', 'weights.h5'),
    'database': database,
}
