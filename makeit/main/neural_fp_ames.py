from __future__ import print_function
from makeit.utils.parsing import input_to_bool, smiles_to_boolean_fps
from makeit.utils.saving import save_model_history, save_model_history_manual
from makeit.utils.parse_cfg import read_config
from makeit.utils.neural_fp import *
from makeit.utils.GraphEmbedding import GraphFP
from keras.models import Sequential, model_from_json
from keras.layers.core import Dense, Dropout, Activation, Masking
from keras.layers.embeddings import Embedding
from keras.layers.recurrent import LSTM
from keras.optimizers import RMSprop, Adam
from keras.utils.visualize_util import plot
import csv
import theano # for looking at computational graph
import rdkit.Chem as Chem
import matplotlib.pyplot as plt
import numpy as np
import datetime
import json
import sys
import os

def build_model(embedding_size = 100, lr = 0.01, optimizer = 'adam', depth = 2, scale_output = 0.05, padding = True):
	'''Generates simple embedding model to use molecular tensor as
	input in order to predict a single-valued output (i.e., yield)'''

	# Base model
	model = Sequential()

	# Add layers
	model.add(GraphFP(embedding_size, sizeAttributeVector() - 1, 
		depth = depth,
		scale_output = scale_output,
		padding = padding,
		activation_inner = 'tanh'))
	print('    model: added GraphFP layer ({} -> {})'.format('mol', embedding_size))
	model.add(Dense(1, activation = 'sigmoid'))
	print('    model: added Dense layer ({} -> {})'.format(embedding_size, 1))

	# Compile
	if optimizer == 'adam':
		optimizer = Adam(lr = lr)
	elif optimizer == 'rmsprop':
		optimizer = RMSProp(lr = lr)
	else:
		print('Can only handle adam or rmsprop optimizers currently')
		quit(1)

	print('compiling...',)
	model.compile(loss = 'mse', optimizer = optimizer)
	print('done')

	return model

def save_model(model, loss, val_loss, fpath = '', config = {}, tstamp = ''):
	'''Saves NN model object and associated information'''

	# Dump data
	with open(fpath + '.json', 'w') as structure_fpath:
		json.dump(model.to_json(), structure_fpath)
	print('...saved structural information')

	# Dump weights
	model.save_weights(fpath + '.h5', overwrite = True)
	print('...saved weights')

	# Dump image
	plot(model, to_file = fpath + '.png')
	print('...saved image')

	# Dump history
	save_model_history_manual(loss, val_loss, fpath + '.hist')
	print ('...saved history')

	# Write to info file
	info_fid = open(fpath + '.info', 'a')
	info_fid.write('{} saved at UTC {}\n\n'.format(fpath, tstamp))
	info_fid.write('Configuration details\n------------\n')
	info_fid.write('  {}\n'.format(config))
	info_fid.close()

	print('...saved model to {}.[json, h5, png, info]'.format(fpath))
	return True

def get_data(data_fpath, training_ratio = 0.9, shuffle_seed = 0):
	'''This is a helper script to read the data .json object and return
	the training and test data sets separately. This is to allow for an
	already-trained model to be evaluated using the test data (i.e., which
	we know it hasn't seen before)'''


	# Load data from json file
	print('reading data...',)
	suppl = Chem.SDMolSupplier(data_fpath)
	print('done')
		
	# Truncate if necessary
	data_indeces = range(len(suppl))
	try:
		truncate_to = int(config['TRAINING']['truncate_to'])
		data_indeces = data_indeces[:truncate_to]
		print('truncated data to first {} samples'.format(truncate_to))
	except:
		pass

	# Get new training_ratio if possible
	try:
		temp = float(config['TRAINING']['ratio'])
		training_ratio = temp
	except:
		pass

	# Get new shuffle seed if possible
	try:
		temp = int(config['TRAINING']['shuffle_seed'])
		shuffle_seed = temp
	except:
		pass

	# Parse data into individual components
	smiles = []
	mols = []
	y = []
	print('processing data...',)
	# Randomize
	np.random.seed(shuffle_seed) # same shuffle every time
	np.random.shuffle(data_indeces)
	for j, i in enumerate(data_indeces):
		if (j % 100) == 99:
			print('  {}/{}'.format(j + 1, len(data_indeces) - 1))
		#try:
		# Molecule first (most likely to fail)
		mol = suppl[i]
		if not mol:
			print('Failed to generate mol for entry {}'.format(i))
			continue
		mols.append(molToGraph(mol).dump_as_tensor())
		nonmutagenic = 'nonmutagen' in suppl.GetItemText(i)
		y.append(1 - int(nonmutagenic)) # 1 if mutagenic
		smiles.append(Chem.MolToSmiles(mol)) # Smiles
		#except:
		#	print('Failed to generate graph for chemical {}'.format(i))

	if int(config['TRAINING']['batch_size']) > 1: # NEED TO PAD
		num_atoms = [x.shape[0] for x in mols]
		max_num_atoms = max(num_atoms)
		print('padding tensors up to N_atoms = {}...'.format(max_num_atoms + 1))
		mols = [padGraphTensor(x, max_num_atoms + 1) for x in mols]
	print('done')

	# Divide up data
	if 'ratio' in config['TRAINING']['data_split']: # split train/notrain
		# Create training/development split
		division = int(len(data) * training_ratio)
		mols_train = mols[:division]
		mols_notrain  = mols[division:]
		y_train = y[:division]
		y_notrain  = y[division:]
		smiles_train = smiles[:division]
		smiles_notrain = smiles[division:]

		# Split notrain up
		mols_val    = mols_notrain[:(len(mols_notrain) / 2)] # first half
		y_val       = y_notrain[:(len(mols_notrain) / 2)] # first half
		smiles_val  = smiles_notrain[:(len(mols_notrain) / 2)] # first half
		mols_test   = mols_notrain[(len(mols_notrain) / 2):] # second half
		y_test      = y_notrain[(len(mols_notrain) / 2):] # second half
		smiles_test = smiles_notrain[(len(mols_notrain) / 2):] # second half

	elif 'cv' in config['TRAINING']['data_split']: # cross-validation
		# Default to first fold of 5-fold cross-validation
		folds = 5
		this_fold = 0
		# Read from config file
		try:
			folds = int(config['TRAINING']['folds'].split('/')[1])
			this_fold = int(config['TRAINING']['folds'].split('/')[0]) - 1
		except:
			pass

		# Get target size of each fold
		N = len(mols)
		target_fold_size = int(np.ceil(float(N) / folds))
		# Split up data
		folded_mols 	= [mols[x:x+target_fold_size]   for x in range(0, N, target_fold_size)]
		folded_y 		= [y[x:x+target_fold_size]      for x in range(0, N, target_fold_size)]
		folded_smiles 	= [smiles[x:x+target_fold_size] for x in range(0, N, target_fold_size)]
		print('Split data into {} folds'.format(folds))
		print('...using fold {}'.format(this_fold + 1))

		# Recombine into validation (this_fold), training (everything else), and testing (same as validation)
		mols_train   = [x for fold in (folded_mols[:this_fold] + folded_mols[(this_fold + 1):])     for x in fold]
		y_train      = [x for fold in (folded_y[:this_fold] + folded_y[(this_fold + 1):])           for x in fold]
		smiles_train = [x for fold in (folded_smiles[:this_fold] + folded_smiles[(this_fold + 1):]) for x in fold]
		# Validation is just this fold
		mols_val     = folded_mols[this_fold]
		y_val        = folded_y[this_fold]
		smiles_val   = folded_smiles[this_fold]
		# Test is a copy of validation
		mols_test    = mols_val[:]
		y_test       = y_val[:]
		smiles_test  = smiles_val[:]

	else:
		print('Must specify a data_split type of "ratio" or "cv"')
		quit(1)

	return (mols_train, y_train, smiles_train,
			mols_val, y_val, smiles_val,
			mols_test, y_test, smiles_test)

def train_model(model, data_fpath = '', nb_epoch = 0, batch_size = 1, lr_func = '0.01', patience = 10):
	'''Trains the model'''

	# Get data from helper function
	(mols_train, y_train, smiles_train,
			mols_val, y_val, smiles_val,
			_, _, _) = get_data(data_fpath)

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
					sloss = model.train_on_batch(single_mol_as_array, single_y_as_array, accuracy = False)
					# print('single loss: {}'.format(sloss))
					this_loss.append(sloss)
				
				# Run through testing set
				print('Testing...')
				for j in range(len(mols_val)):
					single_mol_as_array = np.array(mols_val[j:j+1])
					spred = float(model.predict_on_batch(single_mol_as_array)[0][0][0])
					sloss = (spred - y_val[j]) ** 2
					# print('single loss: {}'.format(sloss))
					this_val_loss.append(sloss)
				
				loss.append(np.mean(this_loss))
				val_loss.append(np.mean(this_val_loss))
				print('loss: {}\tval_loss: {}'.format(loss[i], val_loss[i]))
		
		else: # PADDED VALUES 
			mols = np.vstack((mols_train, mols_val))
			y = np.concatenate((y_train, y_val))
			hist = model.fit(mols, y, 
				nb_epoch = nb_epoch, 
				batch_size = batch_size, 
				validation_split = (1 - float(len(mols_train))/(len(mols_val) + len(mols_train))),
				verbose = 1,
				callbacks = [LearningRateScheduler(lr), EarlyStopping(patience = patience, verbose = 1)])	
			loss = hist.history['loss']
			val_loss = hist.history['val_loss']

	except KeyboardInterrupt:
		print('terminated training early (intentionally)')

	return (model, loss, val_loss)

def test_model(model, data_fpath, fpath, tstamp = '', batch_size = 128):
	'''This function evaluates model performance using test data. The output is
	more meaningful than just the loss function value.'''
	# Create folder to dump testing info to
	try:
		os.makedirs(fpath)
	except: # file exists
		pass
	test_fpath = os.path.join(fpath, tstamp)

	# Get data from helper function
	(mols_train, y_train, smiles_train,
			mols_val, y_val, smiles_val,
			mols_test, y_test, smiles_test) = get_data(data_fpath)

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
		y_train_pred = model.predict(np.array(mols_train), batch_size = batch_size, verbose = 1)[:, 0]
		y_val_pred = model.predict(np.array(mols_val), batch_size = batch_size, verbose = 1)[:, 0]
		y_test_pred = model.predict(np.array(mols_test), batch_size = batch_size, verbose = 1)[:, 0]

	# Save
	with open(test_fpath + '.test', 'w') as fid:
		time_now = datetime.datetime.utcnow()
		fid.write('-- tested at UTC {}\n\n'.format(fpath, time_now))		
		fid.write('test entry\tsmiles\tactual mutagen\tpredicted mutagen\tactual - predicted\n')
		for i in range(len(smiles_test)):
			fid.write('{}\t{}\t{}\t{}\t{}\n'.format(i, 
				smiles_test[i],
				y_test[i], 
				y_test_pred[i],
				y_test[i] - y_test_pred[i]))

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
		Q = (TP + TN) / N
		print('N[{}]\tTP[{}]\tTN[{}]\tFP[{}]\tFN[{}]\tSEN[{}]\tSPEC[{}]\tBAC[{}]\tQ[{}]'.format(
			N, TP, TN, FP, FN, SEN, SPEC, BAC, Q))
		return (N, TP, TN, FP, FN, SEN, SPEC, BAC, Q)

	# Save/report stats
	with open(test_fpath + '.stats', 'w') as fid:
		time_now = datetime.datetime.utcnow()
		fid.write('-- tested at UTC {}\n\n'.format(fpath, time_now))		
		fid.write('\t'.join(['dataset', 'N', 'TP', 'TN', 'FP', 'FN', 'SEN', 'SPEC', 'BAC', 'Q']) + '\n')
		print('Training performance')
		(N, TP, TN, FP, FN, SEN, SPEC, BAC, Q) = score_classification(np.array(y_train) > 0.5, np.array(y_train_pred) > 0.5)
		fid.write('test\t' + '\t'.join(['{}'.format(x) for x in [N, TP, TN, FP, FN, SEN, SPEC, BAC, Q]]))
		print('Validation performance')
		(N, TP, TN, FP, FN, SEN, SPEC, BAC, Q) = score_classification(np.array(y_val) > 0.5, np.array(y_val_pred) > 0.5)
		fid.write('val\t' + '\t'.join(['{}'.format(x) for x in [N, TP, TN, FP, FN, SEN, SPEC, BAC, Q]]))
		print('Testing performance')
		(N, TP, TN, FP, FN, SEN, SPEC, BAC, Q) = score_classification(np.array(y_test) > 0.5, np.array(y_test_pred) > 0.5)
		fid.write('test\t' + '\t'.join(['{}'.format(x) for x in [N, TP, TN, FP, FN, SEN, SPEC, BAC, Q]]))
	
	return

def test_embeddings_demo(model, data_fpath, fpath):
	'''This function tests molecular representations'''
	print('Building images of fingerprint examples')

	# Create folder to dump testing info to
	try:
		os.makedirs(fpath)
	except: # folder exists
		pass
	try:
		fpath = os.path.join(fpath, 'embeddings')
		os.makedirs(fpath)
	except: # folder exists
		pass

	# Get data
	(mols_train, y_train, smiles_train,
			mols_val, y_val, smiles_val,
			mols_test, y_test, smiles_test) = get_data(data_fpath)

	# Define function to test embedding
	tf = K.function([model.layers[0].input], 
		model.layers[0].get_output(train = False))

	# Define function to save image
	def embedding_to_png(embedding, label, fpath):
		fig = plt.figure(figsize=(20,0.5))
		plt.pcolor(embedding, vmin = 0, vmax = 1)
		plt.title('{}'.format(label))
		# cbar = plt.colorbar()
		plt.gca().yaxis.set_visible(False)
		plt.gca().xaxis.set_visible(False)
		plt.xlim([0, 512])
		plt.subplots_adjust(left = 0, right = 1, top = 0.4, bottom = 0)
		plt.savefig(os.path.join(fpath, label) + '.png', bbox_inches = 'tight')
		plt.close(fig)
		plt.clf()
		return

	# Run through training set
	for j in range(len(mols_train)):
		single_mol_as_array = np.array(mols_train[j:j+1])
		embedding = tf([single_mol_as_array])
		embedding_to_png(embedding, smiles_train[j], fpath)
		
		if j == 25:
			break

	# Run through validation set
	for j in range(len(mols_val)):
		single_mol_as_array = np.array(mols_val[j:j+1])

	# Run through testing set
	for j in range(len(mols_test)):
		single_mol_as_array = np.array(mols_test[j:j+1])

	smiles = ''
	while True:
		smiles = raw_input('Enter smiles: ').strip()
		if smiles is 'done':
			break
		try:
			mol = Chem.MolFromSmiles(smiles)
			mol_graph = molToGraph(mol).dump_as_tensor()
			single_mol_as_array = np.array([mol_graph])
			embedding = tf([single_mol_as_array])
			embedding_to_png(embedding, smiles, fpath)
		except:
			print('error saving embedding - was that a SMILES string?')

	return

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

	# Are we loading from an old structure?
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
				padding = int(config['TRAINING']['batch_size']) > 1)
			print('...built untrained model')
		except KeyError:
			print('Missing key in config')
			quit(1)
		# # except ValueError:
		# # 	print('embedding_size must be suitable number')
		# 	quit(1)

	# See if weights exist in this location already
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

	# Testing embeddings?
	try:
		if input_to_bool(config['TESTING']['test_embedding']):
			# Test current embebddings
			test_embeddings_demo(model, config['IO']['data_fpath'])
			quit(1)
	except KeyError:
		pass

	# Train model
	hist = None
	try:
		print('...training model')
		model.optimizer.lr.set_value(float(config['TRAINING']['lr']))
		(model, loss, val_loss) = train_model(model, 
			data_fpath = config['IO']['data_fpath'], 
			nb_epoch = int(config['TRAINING']['nb_epoch']), 
			batch_size = int(config['TRAINING']['batch_size']),
			lr_func = config['TRAINING']['lr'],
			patience = int(config['TRAINING']['patience']))
		print('...trained model')
	except KeyError:
		print('Must specify data_fpath in IO and nb_epoch and batch_size in TRAINING in config')
		quit(1)
#	except ValueError:
#		print('nb_epoch and batch_size must be integers')
#		quit(1)
	except KeyboardInterrupt:
		pass

	# Get the current time
	tstamp = datetime.datetime.utcnow().strftime('%m-%d-%Y_%H-%M')

	# Save for future
	print('...saving model')
	save_model(model, 
		loss,
		val_loss,
		fpath = fpath,
		config = config, 
		tstamp = tstamp)
	print('...saved model')
		
	# Test model
	print('...testing model')
	test_model(model, config['IO']['data_fpath'],
		fpath,
		tstamp = tstamp,
		batch_size = int(config['TRAINING']['batch_size']))
	print('...tested model')
