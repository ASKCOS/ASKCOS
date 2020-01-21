import os
time_zero = 0

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
context_module = neural_network
synth_enumeration = template
retro_enumeration = template
prioritizaton = heuristic
forward_scoring = network

################################################################################
# Define data file locations
################################################################################

data_path = os.path.join(os.path.dirname(__file__),'data')
local_db_dumps = os.path.join(data_path, 'local_db_dumps')
models_path = os.path.join(data_path, 'models')

fingerprint_bits = 256
reaction_fingerprint_bits = 2048


database = 'askcos'

################################################################################
# Define databases (should be nonessential if all local files present)
################################################################################

MONGO_HOST = os.environ.get('MONGO_HOST', 'MONGO_HOST')
MONGO_USER = os.environ.get('MONGO_USER', 'USERNAME')
MONGO_PW = os.environ.get('MONGO_PW', 'PASSWORD')

MONGO = {
    'path': 'mongodb://{}:{}@{}'.format(MONGO_USER, MONGO_PW, MONGO_HOST),
    'id': 27017,
    'connect': False
}

RETRO_TEMPLATES = {
    'file_name': os.path.join(data_path, 'templates', 'retro.templates.json.gz'),
    'database': database,
    'collection': 'retro_templates'
}

FORWARD_TEMPLATES = {
    'file_name': os.path.join(data_path, 'templates', 'forward.templates.json.gz'),
    'database': database,
    'collection': 'forward_templates'
}

REACTIONS = {
    'database': database,
    'collection': 'reactions',
}

CHEMICALS = {
    'file_name': os.path.join(data_path, 'historian', 'chemicals.json.gz'),
    'database': database,
    'collection': 'chemicals',
}

BUYABLES = {
    'file_name': os.path.join(data_path, 'buyables', 'buyables.json.gz'),
    'database': database,
    'collection': 'buyables',
}

SOLVENTS = {
    'file_name': os.path.join(data_path, 'solvents', 'abraham_solvents.pkl'),
    'database': database,
    'collection': 'solvents',
}

# Template-based forward predictor
PREDICTOR = {
    'trained_model_path': os.path.join(models_path, 'template_prioritization', 'forward_scoring'),
    'info': '01-23-17, model trained on 80k Reaxys examples, validated on 10k, tested on 10k. Nh1_200, Nh2_200, Nh3_200, l2_0, Nc_5000, enh_weight_0d1, context_weight_50, opt_adadelta, batch_5, moreFeatures'
}

# Fast filter evaluation
FAST_FILTER_MODEL = {
    'trained_model_path': os.path.join(models_path, 'fast_filter','fast_filter_cleandata.h5'),
}

# Hard coded mincounts to maintain compatibility of the relevance method (weights are numpy matrices)
RELEVANCE_TEMPLATE_PRIORITIZATION = {
    'reaxys': {
        'model_path': os.path.join(models_path, 'template_prioritization', 'reaxys', '1'),
        'output_size': 163723
    }
}

# Different SCScore models that are all functionally similary
SCScore_Prioritiaztion = {
    'trained_model_path_1024bool': os.path.join(models_path, 'scscore', 'model_1024bool.pickle'),
    'trained_model_path_2048bool': os.path.join(models_path, 'scscore', 'model_2048bool.pickle'),
    'trained_model_path_1024uint8': os.path.join(models_path, 'scscore', 'model_1024uint8.pickle')}

MinCost_Prioritiaztion = {
    'trained_model_path': os.path.join(models_path, 'mincost', 'model.hdf5')
}

NEURALNET_CONTEXT_REC = {
    'info_path': os.path.join(models_path, 'context', 'NeuralNet_Cont_Model/'),
    'model_path': os.path.join(models_path, 'context', 'NeuralNet_Cont_Model', 'model.json'),
    'weights_path': os.path.join(models_path, 'context', 'NeuralNet_Cont_Model', 'weights.h5'),
    'database': database,
}

TEMPLATE_FREE_FORWARD_PREDICTOR = {
    'core_model_path': os.path.join(models_path, 'forward_predictor', 'rexgen_direct', 'core_wln_global', 'model-300-3-direct', 'model.ckpt-140000'),
    'rank_model_path': os.path.join(models_path, 'forward_predictor', 'rexgen_direct', 'rank_diff_wln', 'model-core16-500-3-max150-direct-useScores', 'model.ckpt-2400000')
}

SELECTIVITY = {
    'model_path': os.path.join(models_path, 'selectivity', 'model.ckpt-30615')
}
