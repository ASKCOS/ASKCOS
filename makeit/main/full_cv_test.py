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
fpath = 'makeit/models/Abraham_CV_{}'.format(now)
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
	'data_label': 'abraham',
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
	
	# Get data
	data_kwargs['cv_folds'] = '{}/{}'.format(fold_num, N_folds)
	(train, val, test) = get_data_full(**data_kwargs)
	if len(val['mols']) != 0:
		raise ValueError('Validation data should be empty in CV right now!')

	# Fix random subsets for INTERNAL performance testing
	training_indices = range(len(train['mols']))
	internal_trainings = []
	internal_testings = []
	for i in range(3):
		np.random.shuffle(training_indices)
		split = int(len(training_indices) * 0.8)
		internal_trainings.append(
			dict((key, val[:split]) for key, val in train.iteritems())
		) 
		internal_testings.append(
			dict((key, val[split:]) for key, val in train.iteritems())
		)
	print('Generated 3 internal training/testing splits')

	# Randomly select hyperparameters
	param_sets = []
	mse_sets = []
	for i in range(3):
		A = np.power(10., np.random.uniform(-4., -1.))
		decay = np.random.uniform(20., 80.)
		params = {
			'lr_func': 'float({} * np.exp(- epoch / {}))'.format(A, decay),
			'drop_1': np.random.uniform(0., 0.7) ** 2,
			'drop_2': np.random.uniform(0., 0.7) ** 2
		}
		param_sets.append(params)
		mse_sets.append([])
		print('Generated random conditions for trial {}'.format(i))
		print(params)

		# Try for however many internal trainings specified
		for j in range(len(internal_trainings)):
			# Reset model with these parameters
			model = reset_layers.reset(model)
			model.layers[1].p = params['drop_1']
			model.layers[3].p = params['drop_2']

			# Defining training parameters
			train_kwargs = {
				'patience': 1000000,
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

		# Average out MSE from internal tests
		mse_sets[i] = np.mean(mse_sets[i])
		print('Average results of internal trainings, mse = '.format(mse_sets[i]))

	print('Finished all parameter sets, averaged MSEs:')
	print(mse_sets)
	# Find highest ~2~ performing parameter sets
	indices_best = np.array(mse_sets).argsort()[-2:][::-1]

	# Train on full data using best parameters
	residuals = None
	for best, i in enumerate(indices_best):
		params = param_sets[i]

		# Record
		with open(os.path.join(this_fpath, 'model{}.params'.format(best)), 'w') as fid:
			fid.write(str(params))

		# Reset model with these parameters
		model = reset_layers.reset(model)
		model.layers[1].p = params['drop_1']
		model.layers[3].p = params['drop_2']

		# Defining training parameters
		train_kwargs = {
			'patience': 1000000,
			'batch_size': 10,
			'nb_epoch': 100,
			'lr_func': params['lr_func']
		}

		# Train using internal test as validation
		this_data = (train, {}, {})
		(model, loss, val_loss) = train_model(model, this_data, **train_kwargs)

		# Test
		this_data = (train, {}, test)
		(_, _, test_with_resids) = test_model(model, this_data, fpath = this_fpath, 
			tstamp = 'params{}'.format(best), batch_size = 10)

		# Save running average of residuals
		if not residuals:
			residuals = test['residuals']
		else:
			residuals = ((residuals * best) + test['residuals']) / (best + 1)

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