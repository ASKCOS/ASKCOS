from __future__ import print_function
from makeit.utils.parsing import input_to_bool
from makeit.utils.parse_cfg import read_config
from keras.models import model_from_json
import rdkit.Chem as Chem
import matplotlib.pyplot as plt
import datetime
import json
import sys
import os

from makeit.main.core import build_model, train_model, save_model
from makeit.main.test import test_model, test_reactions, test_embeddings_demo
from makeit.main.data import get_data_full

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
		use_old_structure = input_to_bool(config['ARCHITECTURE']['use_existing'])
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
			model = build_model(embedding_size = int(config['ARCHITECTURE']['embedding_size']),
				depth = int(config['ARCHITECTURE']['depth']),
				scale_output = float(config['ARCHITECTURE']['scale_output']),
				padding = int(config['TRAINING']['batch_size']) > 1,
				hidden = int(config['ARCHITECTURE']['hidden']))
			print('...built untrained model')
		except KeyboardInterrupt:
			print('User cancelled model building')
			quit(1)

	###################################################################################
	### LOAD WEIGHTS?
	###################################################################################

	weights_fpath = fpath + '.h5'

	try:
		use_old_weights = input_to_bool(config['TRAINING']['use_existing'])
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
	### CHECK FOR TESTING CONDITIONS
	###################################################################################

	# Testing embeddings?
	try:
		if input_to_bool(config['TESTING']['test_embedding']):
			# Test current embebddings
			test_embeddings_demo(model, get_data, fpath)
			quit(1)
	except KeyError:
		pass

	# Testing activations?
	try:
		if input_to_bool(config['TESTING']['test_activations']):
			# Test current embebddings
			test_activations(model, get_data, fpath, contribution_threshold = 0.5)
			quit(1)
	except KeyError:
		pass

	# Testing reaction embeddings?
	try:
		if input_to_bool(config['TESTING']['test_reactions']):
			test_reactions(model, fpath)
			quit(1)
	except KeyError:
		pass

	###################################################################################
	### DEFINE GET_DATA FUNCTION
	###################################################################################

	def get_data():	return get_data_full(config['IO']['data_label'], config = {})

	###################################################################################
	### TRAIN THE MODEL
	###################################################################################

	# Train model
	try:
		print('...training model')
		(model, loss, val_loss) = train_model(model, 
			get_data = get_data,
			nb_epoch = int(config['TRAINING']['nb_epoch']), 
			batch_size = int(config['TRAINING']['batch_size']),
			lr_func = config['TRAINING']['lr'],
			patience = int(config['TRAINING']['patience']))
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
	test_model(model, get_data,	fpath, tstamp = tstamp,
		batch_size = int(config['TRAINING']['batch_size']))
	print('...tested model')