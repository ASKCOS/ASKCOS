from __future__ import print_function
from qsar.utils.parsing import input_to_bool
from qsar.utils.parse_cfg import read_config
from qsar.utils.neural_fp import sizeAttributeVector
import makeit.utils.reset_layers as reset_layers
from keras.models import model_from_json
import rdkit.Chem as Chem
import matplotlib.pyplot as plt
import datetime
import json
import sys
import os
import time
import numpy as np

from qsar.main.core import build_model
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

	try:
		weights_fpath = config['IO']['weights_fpath']
	except KeyError:
		pass

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

	##############################
	### DEFINE TESTING CONDITIONS
	##############################

	conditions = [[]]
	for i in range(1):#range(sizeAttributeVector() - 1):
		conditions += [
			np.array(i)
		]

	for i, condition in enumerate(conditions):
		print('REMOVING: {}'.format(condition))

		# Get *FULL* dataset
		data_kwargs['training_ratio'] = 1.0
		data_kwargs['cv_folds'] = '1/1'
		data = get_data_full(**data_kwargs)

		# Now filter as needed
		for i in range(3): # for each train, val, test
			for j in range(len(data[i]['mols'])): # for each mol in that list
				data[i]['mols'][j][:, :, condition] = 0.0 # reset that feature index to zero

		###################################################################################
		### TEST MODEL
		###################################################################################
		if type(condition) != type([]):
			stamp = 'reset {}'.format(condition)
		else:
			stamp = 'baseline'
		print('...testing model')
		data_withresiduals = test_model(model, data, fpath, tstamp = stamp,
			batch_size = int(config['TRAINING']['batch_size']))
		print('...tested model')
