import os

# Should we use stereochemistry?
USE_STEREOCHEMISTRY = True

#Output debugging statements
DEBUG = False

#Use highest protocol in pickle
protocol = -1
data_path = os.path.join(os.getcwd(),'data')
fingerprint_bits = 256
pricer_data = os.path.join(data_path,'buyable')
retro_template_data = os.path.join(data_path,'retro-synthetic')
synth_template_data = os.path.join(data_path,'synthetic')
database = 'reaxys_v2'

#Required database names
MONGO = {
        'path': 'mongodb://guest:guest@rmg.mit.edu/admin', 
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
    'collection': 'transforms_retro_v8', # 'lowe' or 'chematica'
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
    'trained_model_path': os.path.join(os.getcwd(),'data/forward_scoring'),
    'info': '01-23-17, model trained on 80k Reaxys examples, validated on 10k, tested on 10k. Nh1_200, Nh2_200, Nh3_200, l2_0, Nc_5000, enh_weight_0d1, context_weight_50, opt_adadelta, batch_5, moreFeatures'
}

CONTEXT_REC = {
    'info_path': os.path.join(data_path,'context', 'RxnID_infoFull.txt'),
    'model_path': os.path.join(data_path,'context', 'fp256noFtr_NN10_BT.pickle'),
    'model_dir': data_path,
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
    'data_loc': os.path.join(data_path, 'context/flow_condition_data4'), 
    'model_loc': os.path.join(data_path, 'context/flow_model4.h5')
    }

FLOW_CONDITIONS4_50 = {
    'database':database,
    'collection' : 'flow_conditions_50',
    'data_loc': os.path.join(data_path, 'context/flow_condition_data4_50'), 
    'model_loc': os.path.join(data_path, 'context/flow_model4_50.h5')
    }