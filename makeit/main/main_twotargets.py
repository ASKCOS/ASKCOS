from __future__ import print_function
from makeit.utils.parsing import input_to_bool
from makeit.utils.parse_cfg import read_config
import makeit.utils.reset_layers as reset_layers
from keras.models import model_from_json
import rdkit.Chem as Chem
import matplotlib.pyplot as plt
import numpy as np
import datetime
import json
import sys
import os
import time

from makeit.main.core import save_model
from makeit.main.data import get_data_full

def build_model(embedding_size = 100, lr = 0.01, optimizer = 'adam', depth = 2, 
	scale_output = 0.05, padding = True, inner_reg_l2 = 0.0, output_reg_l2 = 0.0,
	hidden = 0, loss = 'mse', class_mode = 'categorical', hidden_activation = 'tanh',
	output_activation = 'linear', dr1 = 0.1, dr2 = 0.3):
	'''Generates simple embedding model to use molecular tensor as
	input in order to predict a single-valued output (i.e., yield)

	inputs:
		embedding_size - size of fingerprint for GraphFP layer
		lr - learning rate to use (train_model overwrites this value)
		optimizer - optimization function to use (adam or rmsprop supported)
		depth - depth of the neural fingerprint (i.e., radius)
		scale_output - initial scale for output weights in GraphFP
		padding - whether or not molecular tensors will be padded 
		inner_reg_l2 - l2 regularization parameter for inner GraphFP weights
		output_reg_l2 - l2 regulariztion parameter for outer GraphFP weights 
		hidden - number of hidden tanh nodes after FP (0 is linear)
		loss - loss function as a string (e.g., 'mse')

	outputs:
		model - a Keras model'''

	# Base model
	model = Sequential()

	# Get regularizers
	if inner_reg_l2 > 0:
		inner_regularizer = l2(inner_reg_l2)
	else:
		inner_regularizer = None
	if output_reg_l2 > 0:
		output_regularizer = l2(output_reg_l2)
	else:
		output_regularizer = None

	# Add layers
	model.add(GraphFP(embedding_size, sizeAttributeVector() - 1, 
		depth = depth,
		scale_output = scale_output,
		padding = padding,
		activation_inner = 'tanh',
		inner_regularizer = inner_regularizer,
		output_regularizer = output_regularizer))
	print('    model: added GraphFP layer ({} -> {})'.format('mol', embedding_size))
	if hidden > 0:
		if dr1 > 0:
			model.add(Dropout(dr1))
			print('    model: Added Dropout({})'.format(dr1))
		model.add(Dense(hidden, activation = hidden_activation))
		print('    model: added {} Dense layer (-> {})'.format(hidden_activation, hidden))
		if dr2 > 0:
			model.add(Dropout(dr2))
			print('    model: added Dropout({})'.format(dr2))
	model.add(Dense(2, activation = output_activation))
	print('    model: added lin Dense layer (-> {})'.format(2))

	# Compile
	if optimizer == 'adam':
		optimizer = Adam(lr = lr)
	elif optimizer == 'rmsprop':
		optimizer = RMSprop(lr = lr)
	else:
		print('Can only handle adam or rmsprop optimizers currently')
		quit(1)

	print('compiling...',)
	model.compile(loss = loss, optimizer = optimizer, class_mode = class_mode,
					sample_weight_mode = 'temporal')
	print('done')

	return model

def train_model(model, data, nb_epoch = 0, batch_size = 1, lr_func = '0.01', patience = 10):
	'''Trains the model.

	inputs:
		model - a Keras model
		data - three dictionaries for training,
				validation, and testing separately
		nb_epoch - number of epochs to train for
		batch_size - batch_size to use on the data. This must agree with what was
				specified for data (i.e., if tensors are padded or not)
		lr_func - string which is evaluated with 'epoch' to produce the learning 
				rate at each epoch 
		patience - number of epochs to wait when no progress is being made in 
				the validation loss

	outputs:
		model - a trained Keras model
		loss - list of training losses corresponding to each epoch 
		val_loss - list of validation losses corresponding to each epoch'''

	# Get data from helper function
	(train, val, test) = data
	# Unpack
	mols_train = train['mols']; y_train = train['y']; smiles_train = train['smiles']
	mols_val   = val['mols'];   y_val   = val['y'];   smiles_val   = val['smiles']

	# Create learning rate function
	lr_func_string = 'def lr(epoch):\n    return {}\n'.format(lr_func)
	exec lr_func_string

	# Fit (allows keyboard interrupts in the middle)
	# Because molecular graph tensors are different sizes based on N_atoms, can only do one at a time
	# (alternative is to pad with zeros and try to add some masking feature to GraphFP)
	try:
		loss = []
		val_loss = []

		if batch_size == 1: # DO NOT NEED TO PAD
			wait = 0
			prev_best_val_loss = 99999999
			for i in range(nb_epoch):
				print('Epoch {}/{}, lr = {}'.format(i + 1, nb_epoch, lr(i)))
				this_loss = []
				this_val_loss = []
				model.optimizer.lr.set_value(lr(i))
				
				# Run through training set
				print('Training...')
				training_order = range(len(mols_train))
				np.random.shuffle(training_order)
				for j in training_order:
					single_mol_as_array = np.array(mols_train[j:j+1])
					single_y_as_array = np.array(y_train[j:j+1])
					sample_weights = (~np.isnan(single_y_as_array)).astype(np.float)
					print(single_y_as_array)
					print(sample_weights)
					sloss = model.train_on_batch(single_mol_as_array, single_y_as_array, 
						accuracy = False, sample_weights = sample_weights)
					# print('single loss: {}'.format(sloss))
					this_loss.append(sloss)
				
				# Run through testing set
				print('Testing...')
				for j in range(len(mols_val)):
					single_mol_as_array = np.array(mols_val[j:j+1])
					spred = float(model.predict_on_batch(single_mol_as_array)[0][0])
					if np.isnan(y_val[j][1]):
						sloss = (spred[0] - y_val[j][0]) ** 2
					else:
						sloss = (spred[1] - y_val[j][1]) ** 2
					# print('single loss: {}'.format(sloss))
					this_val_loss.append(sloss)
				
				loss.append(np.mean(this_loss))
				val_loss.append(np.mean(this_val_loss))
				print('mse loss: {}\tmse val_loss: {}'.format(loss[i], val_loss[i]))

				# Check progress
				if np.mean(this_val_loss) < prev_best_val_loss:
					wait = 0
					prev_best_val_loss = np.mean(this_val_loss)
					if patience == -1:
						model.save_weights('best.h5', overwrite=True)
				else:
					wait = wait + 1
					print('{} epochs without val_loss progress'.format(wait))
					if wait == patience:
						print('stopping early!')
						break
			if patience == -1:
				model.load_weights('best.h5')

		else: # PADDED VALUES 
			callbacks = [LearningRateScheduler(lr)]
			if patience != -1:
				callbacks.append(EarlyStopping(patience = patience, verbose = 1))

			if mols_val:
				mols = np.vstack((mols_train, mols_val))
				y = np.concatenate((y_train, y_val))
				sample_weights = ~np.isnan(y)
				hist = model.fit(mols, y, 
					nb_epoch = nb_epoch, 
					batch_size = batch_size, 
					validation_split = (1 - float(len(mols_train))/(len(mols_val) + len(mols_train))),
					verbose = 1,
					callbacks = callbacks,
					sample_weights = sample_weights)	
			else:
				sample_weights = ~np.isnan(np.array(y_train))
				hist = model.fit(np.array(mols_train), np.array(y_train), 
					nb_epoch = nb_epoch, 
					batch_size = batch_size, 
					verbose = 1,
					callbacks = callbacks,
					sample_weights = sample_weights)	
			
			loss = []; val_loss = []
			if 'loss' in hist.history: loss = hist.history['loss']
			if 'val_loss' in hist.history: val_loss = hist.history['val_loss']

	except KeyboardInterrupt:
		print('terminated training early (intentionally)')

	return (model, loss, val_loss)

def test_model(model, data, fpath, tstamp = 'no_time', batch_size = 128):
	'''This function evaluates model performance using test data. The output is
	more meaningful than just the loss function value.

	inputs:
		model - the trained Keras model
		data - three dictionaries for training,
					validation, and testing data. Each dictionary should have
					keys of 'mol', a molecular tensor, 'y', the target output, 
					and 'smiles', the SMILES string of that molecule
		fpath - folderpath to save test data to, will be appended with '/tstamp.test'
		tstamp - timestamp to add to the testing
		batch_size - batch_size to use while testing'''
	
	# Create folder to dump testing info to
	try:
		os.makedirs(fpath)
	except: # file exists
		pass
	test_fpath = os.path.join(fpath, tstamp)

	# Unpack
	(train, val, test) = data
	# Unpack
	mols_train = train['mols']; y_train = train['y']; smiles_train = train['smiles']
	mols_val   = val['mols'];   y_val   = val['y'];   smiles_val   = val['smiles']
	mols_test  = test['mols'];  y_test  = test['y'];  smiles_test  = test['smiles']
	y_label = test['y_label']

	# Fit (allows keyboard interrupts in the middle)
	# Because molecular graph tensors are different sizes based on N_atoms, can only do one at a time
	# (alternative is to pad with zeros and try to add some masking feature to GraphFP)
	y_train_pred = []
	y_val_pred = []
	y_test_pred = []

	if batch_size == 1: # UNEVEN TENSORS, ONE AT A TIME PREDICTION
		# Run through training set
		for j in range(len(mols_train)):
			single_mol_as_array = np.array(mols_train[j:j+1])
			single_y_as_array = np.array(y_train[j:j+1])
			spred = float(model.predict_on_batch(single_mol_as_array)[0][0][0])
			y_train_pred.append(spred)

		# Run through validation set
		for j in range(len(mols_val)):
			single_mol_as_array = np.array(mols_val[j:j+1])
			spred = float(model.predict_on_batch(single_mol_as_array)[0][0][0])
			y_val_pred.append(spred)

		# Run through testing set
		for j in range(len(mols_test)):
			single_mol_as_array = np.array(mols_test[j:j+1])
			spred = float(model.predict_on_batch(single_mol_as_array)[0][0][0])
			y_test_pred.append(spred)

	else: # PADDED
		y_train_pred = np.array([]); y_val_pred = np.array([]); y_test_pred = np.array([])
		if mols_train: y_train_pred = model.predict(np.array(mols_train), batch_size = batch_size, verbose = 1)[:, 0]
		if mols_val: y_val_pred = model.predict(np.array(mols_val), batch_size = batch_size, verbose = 1)[:, 0]
		if mols_test: y_test_pred = model.predict(np.array(mols_test), batch_size = batch_size, verbose = 1)[:, 0]

	def round3(x):
		return int(x * 1000) / 1000.0

	def parity_plot(true, pred, set_label):
		if len(true) == 0:
			print('skipping parity plot for empty dataset')
			return

		# Calculate some stats
		min_y = np.min((true, pred))
		max_y = np.max((true, pred))
		mse = stats.mse(true, pred)
		mae = stats.mae(true, pred)
		q = stats.q(true, pred)
		(r2, a) = stats.linreg(true, pred) # predicted v observed
		(r2p, ap) = stats.linreg(pred, true) # observed v predicted

		# Create parity plot
		plt.scatter(true, pred, alpha = 0.5)
		plt.xlabel('Actual {}'.format(y_label))
		plt.ylabel('Predicted {}'.format(y_label))
		plt.title('Parity plot for {} ({} set, N = {})'.format(y_label, set_label, len(true)) + 
			'\nMSE = {}, MAE = {}, q = {}'.format(round3(mse), round3(mae), round3(q)) + 
			'\na = {}, r^2 = {}'.format(round3(a), round3(r2)) + 
			'\na` = {}, r^2` = {}'.format(round3(ap), round3(r2p)))
		plt.grid(True)
		plt.plot(true, true * a, 'r--')
		plt.axis([min_y, max_y, min_y, max_y])	
		plt.savefig(test_fpath + ' {}.png'.format(set_label), bbox_inches = 'tight')
		plt.clf()

		# Print
		print('{}:'.format(set_label))
		print('  mse = {}, mae = {}'.format(mse, mae))
		print('  q = {}'.format(q))
		print('  r2 through origin = {} (pred v. true), {} (true v. pred)'.format(r2, r2p))
		print('  slope through origin = {} (pred v. true), {} (true v. pred)'.format(a[0], ap[0]))

	def score_classification(actual, predicted):
		'''Given two numpy boolean vectors, calculates various performance measures'''
		if actual.shape != predicted.shape:
			print('Shapes of inputs do not match in score_classification')
			return

		TP = 0
		TN = 0
		FP = 0
		FN = 0
		N = len(actual)
		for i in range(N):
			if actual[i]:
				if predicted[i]:
					TP += 1
				else:
					FN += 1
			else:
				if predicted[i]:
					FP += 1
				else:
					TN += 1
		SEN = float(TP) / (TP + FN)
		SPEC = float(TN) / (TN + FP)
		BAC = float(SEN + SPEC) / 2.0
		Q = (TP + TN) / float(N)
		print('N[{}]\tTP[{}]\tTN[{}]\tFP[{}]\tFN[{}]\tSEN[{}]\tSPEC[{}]\tBAC[{}]\tQ[{}]'.format(
			N, TP, TN, FP, FN, SEN, SPEC, BAC, Q))
		return (N, TP, TN, FP, FN, SEN, SPEC, BAC, Q)

	# Save
	with open(test_fpath + '.test', 'w') as fid:
		fid.write('{} tested at UTC {}, predicting {}\n\n'.format(fpath, tstamp, y_label))		
		fid.write('test entry\tsmiles\tactual\tpredicted\tactual - predicted\n')
		for i in range(len(smiles_test)):
			fid.write('{}\t{}\t{}\t{}\t{}\n'.format(i, 
				smiles_test[i],
				y_test[i], 
				y_test_pred[i],
				y_test[i] - y_test_pred[i]))

	# Is this a boolean prediction?
	if len(np.unique(y_test)) == 2:
		
		# Save/report stats
		with open(test_fpath + '.stats', 'w') as fid:
			time_now = datetime.datetime.utcnow()
			fid.write('-- tested at UTC {}\n\n'.format(fpath, time_now))		
			fid.write('\t'.join(['dataset', 'N', 'TP', 'TN', 'FP', 'FN', 'SEN', 'SPEC', 'BAC', 'Q']) + '\n')
			print('Training performance')
			(N, TP, TN, FP, FN, SEN, SPEC, BAC, Q) = score_classification(np.array(y_train) > 0.5, np.array(y_train_pred) > 0.5)
			fid.write('test\t' + '\t'.join(['{}'.format(x) for x in [N, TP, TN, FP, FN, SEN, SPEC, BAC, Q]]) + '\n')
			print('Validation performance')
			(N, TP, TN, FP, FN, SEN, SPEC, BAC, Q) = score_classification(np.array(y_val) > 0.5, np.array(y_val_pred) > 0.5)
			fid.write('val\t' + '\t'.join(['{}'.format(x) for x in [N, TP, TN, FP, FN, SEN, SPEC, BAC, Q]]) + '\n')
			print('Testing performance')
			(N, TP, TN, FP, FN, SEN, SPEC, BAC, Q) = score_classification(np.array(y_test) > 0.5, np.array(y_test_pred) > 0.5)
			fid.write('test\t' + '\t'.join(['{}'.format(x) for x in [N, TP, TN, FP, FN, SEN, SPEC, BAC, Q]]) + '\n')

	else:

		# Create plots for datasets
		if y_train: 
			try: parity_plot(y_train, y_train_pred, 'train')
			except: pass
		if y_val: 
			try: parity_plot(y_val, y_val_pred, 'val')
			except: pass
		if y_test: 
			try: parity_plot(y_test, y_test_pred, 'test')
			except: pass

	train['residuals'] = np.array(y_train) - np.array(y_train_pred)
	val['residuals'] = np.array(y_val) - np.array(y_val_pred)
	test['residuals'] = np.array(y_test) - np.array(y_test_pred)

	return (train, val, test)

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
		
		# GET FIRST SET OF DATA
		data_kwargs['data_label'] = 'abraham'
		data_oct = get_data_full(**data_kwargs)
		# GET SECOND SET OF DATA
		data_kwargs['data_label'] = 'delaney'
		data_aq = get_data_full(**data_kwargs)

		# Need to merge datasets
		(train_oct, val_oct, test_oct) = data_oct
		(train_aq, val_aq, test_aq) = data_aq
		train = {}; val = {}; test = {}
		train['mols'] = train_oct['mols'] + train_aq['mols']
		train['smiles'] = train_oct['smiles'] + train_aq['smiles']
		train['y'] = [[x,np.nan] for x in train_oct['y']] + 
		             [[np.nan,x] for x in train_aq['y']]
		print(len(train['mols']))
		print(len(train['y']))
		val['mols'] = val_oct['mols'] + val_aq['mols']
		val['smiles'] = val_oct['smiles'] + val_aq['smiles']
		val['y'] = [[x,np.nan] for x in val_oct['y']] + 
		             [[np.nan,x] for x in val_aq['y']]
		print(len(val['mols']))
		print(len(val['y']))
		test['mols'] = test_oct['mols'] + test_aq['mols']
		test['smiles'] = test_oct['smiles'] + test_aq['smiles']
		test['y'] = [[x,np.nan] for x in test_oct['y']] + 
		             [[np.nan,x] for x in test_aq['y']]
		print(len(test['mols']))
		print(len(test['y']))

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

		#######################
		### RESET MODEL WEIGHTS
		#######################

		model = reset_layers.reset(model)