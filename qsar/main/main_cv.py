from __future__ import print_function
from qsar.utils.parsing import input_to_bool
from qsar.utils.parse_cfg import read_config
import makeit.utils.reset_layers as reset_layers
from keras.models import model_from_json
import rdkit.Chem as Chem
import matplotlib.pyplot as plt
import datetime
import json
import sys
import os
import time

from qsar.main.core import build_model, train_model, save_model
from qsar.main.test import test_model, test_activations, test_embeddings_demo, test_predictions
from qsar.main.data import get_data_full

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print('Usage: {} "settings.cfg"'.format(sys.argv[0]))
		quit(1)

	# Load settings
	try:
		config = read_config(sys.argv[1])
	except:
		print('Could not read config file {}'.format(sys.argv[1]))
		quit(1)

	# Get model label
	try:
		fpath = config['IO']['model_fpath']
	except KeyError:
		print('Must specify model_fpath in IO in config')
		quit(1)

	###################################################################################
	### LOAD STRUCTURE OR BUILD MODEL
	###################################################################################

	structure_fpath = fpath + '.json'
	try:
		use_old_structure = input_to_bool(config['IO']['use_existing_model'])
	except KeyError:
		print('Must specify whether or not to use existing model architecture')
		quit(1)
	if use_old_structure and os.path.isfile(structure_fpath):
		# Load model
		with open(structure_fpath, 'r') as structure_fid:
			print('...loading model architecture')
			model = model_from_json(json.load(structure_fid))
			print('...loaded structural information')
	elif use_old_structure and not os.path.isfile(structure_fpath):
		print('Model not found at specified path {}'.format(structure_fpath))
		quit(1)
	else:
		# Build model
		print('...building model')
		try:
			kwargs = config['ARCHITECTURE']
			del kwargs['__name__'] #  from configparser
			if 'batch_size' in config['TRAINING']:
				kwargs['padding'] = int(config['TRAINING']['batch_size']) > 1
			if 'embedding_size' in kwargs: 
				kwargs['embedding_size'] = int(kwargs['embedding_size'])
			if 'hidden' in kwargs: 
				kwargs['hidden'] = int(kwargs['hidden'])
			if 'depth' in kwargs: 
				kwargs['depth'] = int(kwargs['depth'])
			if 'scale_output' in kwargs: 
				kwargs['scale_output'] = float(kwargs['scale_output'])
			if 'dr1' in kwargs:
				kwargs['dr1'] = float(kwargs['dr1'])
			if 'dr2' in kwargs:
				kwargs['dr2'] = float(kwargs['dr2'])
			if 'output_size' in kwargs:
				kwargs['output_size'] = int(kwargs['output_size'])
				
			model = build_model(**kwargs)
			print('...built untrained model')
		except KeyboardInterrupt:
			print('User cancelled model building')
			quit(1)

	###################################################################################
	### LOAD WEIGHTS?
	###################################################################################

	weights_fpath = fpath + '.h5'

	try:
		use_old_weights = input_to_bool(config['IO']['use_existing_weights'])
	except KeyError:
		print('Must specify whether or not to use existing model weights')
		quit(1)

	if use_old_weights and os.path.isfile(weights_fpath):
		# Load weights
		model.load_weights(weights_fpath)
		print('...loaded weight information')
	elif use_old_weights and not os.path.isfile(weights_fpath):
		print('Weights not found at specified path {}'.format(weights_fpath))
		quit(1)
	else:
		# New weights will be used anyway
		pass
	
	###################################################################################
	### DEFINE DATA 
	###################################################################################

	data_kwargs = config['DATA']
	if '__name__' in data_kwargs:
		del data_kwargs['__name__'] #  from configparser
	if 'batch_size' in config['TRAINING']:
		data_kwargs['batch_size'] = int(config['TRAINING']['batch_size'])
	if 'shuffle_seed' in data_kwargs:
		data_kwargs['shuffle_seed'] = int(data_kwargs['shuffle_seed'])
	else:
		data_kwargs['shuffle_seed'] = int(time.time())
	if 'truncate_to' in data_kwargs:
		data_kwargs['truncate_to'] = int(data_kwargs['truncate_to'])
	if 'training_ratio' in data_kwargs:
		data_kwargs['training_ratio'] = float(data_kwargs['training_ratio'])

	if 'cv_folds' in data_kwargs:
		if '<this_fold>' in data_kwargs['cv_folds']:
			cv_folds = data_kwargs['cv_folds']
			total_folds = int(cv_folds.split('/')[1])
			all_cv_folds = ['{}/{}'.format(i + 1, total_folds) for i in range(total_folds)]
		else:
			all_cv_folds = [data_kwargs['cv_folds']]

	# Iterate through all
	ref_fpath = fpath
	for cv_fold in all_cv_folds:
		print('Using CV fold {}'.format(cv_fold))
		data_kwargs['cv_folds'] = cv_fold
		fpath = ref_fpath.replace('<this_fold>', cv_fold.split('/')[0])
		data = get_data_full(**data_kwargs)

		###################################################################################
		### CHECK FOR TESTING CONDITIONS
		###################################################################################

		# Testing embeddings?
		try:
			if input_to_bool(config['TESTING']['test_embedding']):
				# Test current embebddings
				test_embeddings_demo(model, data, fpath)
				quit(1)
		except KeyError:
			pass

		# Testing activations?
		try:
			if input_to_bool(config['TESTING']['test_activations']):
				# Test current embebddings
				test_activations(model, data, fpath, contribution_threshold = 0.5)
				quit(1)
		except KeyError:
			pass

		###################################################################################
		### TRAIN THE MODEL
		###################################################################################

		# Train model
		try:
			print('...training model')
			kwargs = config['TRAINING']
			if '__name__' in kwargs:
				del kwargs['__name__'] #  from configparser
			if 'nb_epoch' in kwargs:
				kwargs['nb_epoch'] = int(kwargs['nb_epoch'])
			if 'batch_size' in kwargs:
				kwargs['batch_size'] = int(kwargs['batch_size'])
			if 'patience' in kwargs:
				kwargs['patience'] = int(kwargs['patience'])
			(model, loss, val_loss) = train_model(model, data, **kwargs)
			print('...trained model')
		except KeyboardInterrupt:
			pass

		###################################################################################
		### SAVE MODEL
		###################################################################################

		# Get the current time
		tstamp = datetime.datetime.utcnow().strftime('%m-%d-%Y_%H-%M')
		print('...saving model')
		save_model(model, 
			loss,
			val_loss,
			fpath = fpath,
			config = config, 
			tstamp = tstamp)
		print('...saved model')

		###################################################################################
		### TEST MODEL
		###################################################################################

		print('...testing model')
		data_withresiduals = test_model(model, data, fpath, tstamp = tstamp,
			batch_size = int(config['TRAINING']['batch_size']))
		print('...tested model')

		###################################################################################
		### TEST PREDICTIONS
		###################################################################################

		try:
			if input_to_bool(config['TESTING']['test_predictions']):
				test_predictions(model, data_withresiduals, fpath, tstamp = tstamp)
		except KeyError:
			pass

		#######################
		### RESET MODEL WEIGHTS
		#######################

		model = reset_layers.reset(model)