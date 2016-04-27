from __future__ import print_function
from makeit.utils.parsing import input_to_bool
from makeit.utils.parse_cfg import read_config
import makeit.utils.reset_layers as reset_layers
from keras.models import model_from_json
import rdkit.Chem as Chem
import matplotlib.pyplot as plt
import datetime
import time
import json
import sys
import os
import numpy as np

from makeit.main.core import build_model, train_model, save_model
from makeit.main.test import test_model, test_reactions, test_activations, test_embeddings_demo, test_predictions
from makeit.main.data import get_data_full

# Saving to...
now = int(time.time())
fpath = 'makeit/models/Delaney_CV_{}'.format(now)
os.makedirs(fpath)

# Define fixed model architecture and build
model_kwargs = {
	'padding': True,   # if batch_size > 1
	'hidden': 50,      # number of hidden nodes
	'depth': 5,        # radius for blending
}
model = build_model(**model_kwargs)
print('...built untrained model')

# Define data
seed = now # use time in seconds for random seed
data_kwargs = {
	'data_label': 'delaney',
	'batch_size': 10,
	'shuffle_seed': seed,
	'data_split': 'cv',
	'training_ratio': 1.0
}

# Go through folds now
N_folds = 5
test_mses = np.zeros((N_folds,))
for fold_num in range(1, N_folds + 1):
	print('On fold {}'.format(fold_num))
	this_fpath = os.path.join(fpath, 'fold{}'.format(fold_num))
	os.makedirs(this_fpath)
	temp_fpath = os.path.join(this_fpath, 'temp')
	os.makedirs(temp_fpath)
	
	# Get data
	data_kwargs['cv_folds'] = '{}/{}'.format(fold_num, N_folds)
	(train, val, test) = get_data_full(**data_kwargs)
	if len(val['mols']) != 0:
		raise ValueError('Validation data should be empty in CV right now!')

	# Fix random subsets for INTERNAL performance testing
	training_indices = range(len(train['mols']))
	internal_trainings = []
	internal_testings = []
	N_internals = 3
	for i in range(N_internals):
		np.random.shuffle(training_indices)
		split = int(len(training_indices) * 0.8)
		internal_trainings.append(
			dict((key, val[:split]) for key, val in train.iteritems())
		) 
		internal_testings.append(
			dict((key, val[split:]) for key, val in train.iteritems())
		)
	print('Generated {} internal training/testing splits'.format(N_internals))

	# Randomly select hyperparameters
	param_sets = []
	mse_sets = []
	N_trials = 10
	for i in range(N_trials):
		A = np.power(10., np.random.uniform(-2.5, -1.5))
		decay = np.random.uniform(20., 50.)
		params = {
			'lr_func': 'float({} * np.exp(- epoch / {}))'.format(A, decay),
			'drop_1': 0.0, #np.random.uniform(0., 0.01) ** 3,
			'drop_2': np.random.uniform(0., 0.5),
		}
		param_sets.append(params)
		mse_sets.append([])
		print('Generated random conditions for trial {}/{}'.format(i+1, N_trials))
		print(params)

		# Try for however many internal trainings specified
		for j in range(len(internal_trainings)):
			# Reset model with these parameters
			model = reset_layers.reset(model)
			model.layers[1].p = params['drop_1']
			model.layers[3].p = params['drop_2']

			# Defining training parameters
			train_kwargs = {
				'patience': -1,
				'batch_size': 10,
				'nb_epoch': 150,
				'lr_func': params['lr_func']
			}

			# Train using internal test as validation
			this_data = (internal_trainings[j], internal_testings[j], {})
			(model, loss, val_loss) = train_model(model, this_data, **train_kwargs)

			# Save MSE
			if not val_loss: val_loss = [99999999]
			mse_sets[i].append(val_loss[-1])
			print('  internal training {}/{} -> mse = {}'.format(j+1, len(internal_trainings), val_loss[-1]))

			# Save weights to params###_model###.h5
			model.save_weights(os.path.join(temp_fpath, 'params{}_model{}.h5'.format(i, j)))
			print('Saved weights at {}'.format(os.path.join(temp_fpath, 'params{}_model{}.h5'.format(i, j))))

		# Average out MSE from internal tests
		mse_sets[i] = np.mean(mse_sets[i])
		print('Average results of internal trainings, mse = {}'.format(mse_sets[i]))

	print('Finished all parameter sets, averaged MSEs:')
	print(mse_sets)

	# Find highest performing parameter sets (LOWEST error)
	N_best = 2
	indices_best = np.array(mse_sets).argsort()[:N_best]

	# Record
	with open(os.path.join(this_fpath, 'param_performance.txt'), 'w') as fid:
		for i in range(len(mse_sets)):
			using = '*** using in final prediction' if i in indices_best else ''
			fid.write('{}\t{}\t{}\n'.format(mse_sets[i], param_sets[i], using))

	# Train on full data using best parameters
	residuals = None
	num_avg = 0
	for best, i in enumerate(indices_best):
		params = param_sets[i]

		# Record
		with open(os.path.join(this_fpath, 'params{}.params'.format(best)), 'w') as fid:
			fid.write(str(params))

		# Reset model with these parameters
		model = reset_layers.reset(model)
		model.layers[1].p = params['drop_1']
		model.layers[3].p = params['drop_2']

		# # Defining training parameters
		# train_kwargs = {
		# 	'patience': -1,
		# 	'batch_size': 10,
		# 	'nb_epoch': 150,
		# 	'lr_func': params['lr_func']
		# }

		# # Train using internal test as validation
		# a = np.array([])
		# dummy_data = {'mols':a, 'y':a, 'smiles':a}
		# this_data = (train, dummy_data, dummy_data)
		# (model, loss, val_loss) = train_model(model, this_data, **train_kwargs)

		# For each of our internally-trained models:
		for j in range(len(internal_trainings)):
			# Load saved weights
			model.load_weights(os.path.join(temp_fpath, 'params{}_model{}.h5'.format(i, j)))

			# Test
			a = np.array([])
			dummy_data = {'mols':a, 'y':a, 'smiles':a}
			this_data = (dummy_data, dummy_data, test)
			(_, _, test_with_resids) = test_model(model, this_data, fpath = temp_fpath, 
				tstamp = 'params{}_model{}'.format(i, j), batch_size = train_kwargs['batch_size'])

			# Save running average of residuals
			if type(residuals) is type(None):
				residuals = test['residuals']
			else:
				residuals = ((residuals * num_avg) + test['residuals']) / (num_avg + 1)
			num_avg += 1

	# Use averaged residuals to back out averaged test file
	with open(os.path.join(this_fpath, 'averaged.test'), 'w') as fid:
		fid.write('test entry\tsmiles\tactual\tpredicted\tactual - predicted\n')
		for i in range(len(test['smiles'])):
			fid.write('{}\t{}\t{}\t{}\t{}\n'.format(i, 
				test['smiles'][i],
				test['y'][i], 
				test['y'][i] - residuals[i],
				residuals[i]))

	# Save
	test_mses[fold_num - 1] = np.mean(residuals ** 2)
	print('Overall fold MSE: {}'.format(np.mean(residuals ** 2)))

# Report overall MSE
print('Overall MSE: {}'.format(np.mean(test_mses)))
