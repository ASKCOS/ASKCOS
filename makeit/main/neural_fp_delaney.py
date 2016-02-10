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

def build_model(embedding_size = 100, lr = 0.01, optimizer = 'adam', depth = 2, 
	scale_output = 0.05, padding = True):
	'''Generates simple embedding model to use reaction smiles as
	input in order to predict a single-valued output (i.e., yield)'''

	# Base model
	model = Sequential()

	# Add layers
	model.add(GraphFP(embedding_size, sizeAttributeVector() - 1, 
		depth = depth,
		scale_output = scale_output,
		padding = padding))
	print('    model: added GraphFP layer ({} -> {})'.format('mol', embedding_size))
	model.add(Dense(1))
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

def get_data(data_fpath, training_ratio = 0.9):
	'''This is a helper script to read the data .json object and return
	the training and test data sets separately. This is to allow for an
	already-trained model to be evaluated using the test data (i.e., which
	we know it hasn't seen before)'''


	# Load data from json file
	print('reading data...',)
	data = []
	with open(data_fpath, 'r') as data_fid:
		reader = csv.reader(data_fid, delimiter = ',', quotechar = '"')
		for row in reader:
			data.append(row)
	print('done')
		
	# Truncate if necessary
	try:
		truncate_to = int(config['TRAINING']['truncate_to'])
		data = data[:truncate_to]
		print('truncated data to first {} samples'.format(truncate_to))
	except:
		pass

	# Get new training_ratio if possible
	try:
		temp = float(config['TRAINING']['ratio'])
		training_ratio = temp
	except:
		pass

	# Parse data into individual components
	smiles = []
	mols = []
	y = []
	print('processing data...',)
	# Randomize
	np.random.seed(0)
	np.random.shuffle(data)
	for i, row in enumerate(data):
		if i == 0:
			continue
		elif (i % 100) == 99:
			print('  {}/{}'.format(i + 1, len(data) - 1))
		try:
			# Molecule first (most likely to fail)
			mols.append(molToGraph(Chem.MolFromSmiles(row[3])).dump_as_tensor())
			y.append(float(row[1])) # Measured log(solubility M/L)
			smiles.append(row[3]) # Smiles
		except:
			print('Failed to generate graph for {}'.format(row[3]))

	if int(config['TRAINING']['batch_size']) > 1: # NEED TO PAD
		num_atoms = [x.shape[0] for x in mols]
		max_num_atoms = max(num_atoms)
		print('padding tensors up to N_atoms = {}...'.format(max_num_atoms + 1))
		mols = [padGraphTensor(x, max_num_atoms + 1) for x in mols]
	print('done')

	# Create training/development split
	division = int(len(data) * training_ratio)
	mols_train = mols[:division]
	mols_notrain  = mols[division:]
	y_train = y[:division]
	y_notrain  = y[division:]
	smiles_train = smiles[:division]
	smiles_notrain = smiles[division:]

	return (mols_train, y_train, mols_notrain, y_notrain, training_ratio, smiles_train, smiles_notrain)

def train_model(model, data_fpath = '', nb_epoch = 0, batch_size = 1):
	'''Trains the model'''

	# Get data from helper function
	data = get_data(data_fpath)
	mols_train = data[0]
	y_train = data[1]
	mols_notrain = data[2]
	y_notrain = data[3]
	training_ratio = data[4]
	
	# Split notrain up
	mols_val    = mols_notrain[:(len(mols_notrain) / 2)] # first half
	y_val       = y_notrain[:(len(mols_notrain) / 2)] # first half

	# Fit (allows keyboard interrupts in the middle)
	# Because molecular graph tensors are different sizes based on N_atoms, can only do one at a time
	# (alternative is to pad with zeros and try to add some masking feature to GraphFP)
	try:
		loss = []
		val_loss = []

		if batch_size == 1: # DO NOT NEED TO PAD
			for i in range(nb_epoch):
				print('Epoch {}/{}'.format(i + 1, nb_epoch))
				this_loss = []
				this_val_loss = []

				# Shuffle training data

				
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
				verbose = 1)	
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
	data = get_data(data_fpath)
	mols_train = data[0]
	y_train = data[1]
	mols_notrain = data[2]
	y_notrain = data[3]
	training_ratio = data[4]
	smiles_train = data[5]
	smiles_notrain = data[6]

	# Split notrain up
	mols_val    = mols_notrain[:(len(mols_notrain) / 2)] # first half
	y_val       = y_notrain[:(len(mols_notrain) / 2)] # first half
	smiles_val  = smiles_notrain[:(len(mols_notrain) / 2)] # first half
	mols_test   = mols_notrain[(len(mols_notrain) / 2):] # second half
	y_test      = y_notrain[(len(mols_notrain) / 2):] # second half
	smiles_test = smiles_notrain[(len(mols_notrain) / 2):] # second half

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
		fid.write('test entry\tsmiles\tactual log(sol)\tpredicted log(sol)\tactual - predicted\n')
		for i in range(len(smiles_test)):
			fid.write('{}\t{}\t{}\t{}\t{}\n'.format(i, 
				smiles_test[i],
				y_test[i], 
				y_test_pred[i],
				y_test[i] - y_test_pred[i]))

	# Create parity plot
	plt.scatter(y_train, y_train_pred, alpha = 0.5)
	plt.xlabel('Actual log(solubility (M))')
	plt.ylabel('Predicted log(solubility (M))')
	plt.title('Parity plot for solubility prediction (train set)\nMSE({}) = {}'.format(
		len(y_train), np.mean((np.array(y_train) - np.array(y_train_pred)) ** 2)))
	plt.grid(True)
	min_y = min(plt.axis())
	max_y = max(plt.axis())
	plt.axis([min_y, max_y, min_y, max_y])	
	plt.savefig(test_fpath + ' train.png', bbox_inches = 'tight')
	plt.clf()

	# Create parity plot
	plt.scatter(y_val, y_val_pred, alpha = 0.5)
	plt.xlabel('Actual log(solubility (M))')
	plt.ylabel('Predicted log(solubility (M))')
	plt.title('Parity plot for solubility prediction (validation set)\nMSE({}) = {}'.format(
		len(y_val), np.mean((np.array(y_val) - np.array(y_val_pred)) ** 2)))
	plt.grid(True)
	min_y = min(plt.axis())
	max_y = max(plt.axis())
	plt.axis([min_y, max_y, min_y, max_y])
	plt.savefig(test_fpath + ' val.png', bbox_inches = 'tight')
	plt.clf()

	# Create parity plot
	plt.scatter(y_test, y_test_pred, alpha = 0.5)
	plt.xlabel('Actual log(solubility (M))')
	plt.ylabel('Predicted log(solubility (M))')
	plt.title('Parity plot for solubility prediction (test set)\nMSE({}) = {}'.format(
		len(y_test), np.mean((np.array(y_test) - np.array(y_test_pred)) ** 2)))
	plt.grid(True)
	min_y = min(plt.axis())
	max_y = max(plt.axis())
	plt.axis([min_y, max_y, min_y, max_y])	
	plt.savefig(test_fpath + ' test.png', bbox_inches = 'tight')
	plt.clf()

	return

def test_embeddings_demo(model, data_fpath, fpath):
	'''This function tests molecular representations'''
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
	data = get_data(data_fpath)
	mols_train = data[0]
	y_train = data[1]
	mols_notrain = data[2]
	y_notrain = data[3]
	training_ratio = data[4]
	smiles_train = data[5]
	smiles_notrain = data[6]

	# Split notrain up
	mols_val    = mols_notrain[:(len(mols_notrain) / 2)] # first half
	y_val       = y_notrain[:(len(mols_notrain) / 2)] # first half
	smiles_val  = smiles_notrain[:(len(mols_notrain) / 2)] # first half
	mols_test   = mols_notrain[(len(mols_notrain) / 2):] # second half
	y_test      = y_notrain[(len(mols_notrain) / 2):] # second half
	smiles_test = smiles_notrain[(len(mols_notrain) / 2):] # second half


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

	return


def test_activations(model, data_fpath, fpath):
	'''This function tests activation of the output'''
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
	data = get_data(data_fpath)
	mols_train = data[0]
	y_train = data[1]
	mols_notrain = data[2]
	y_notrain = data[3]
	training_ratio = data[4]
	smiles_train = data[5]
	smiles_notrain = data[6]

	# Split notrain up
	mols_val    = mols_notrain[:(len(mols_notrain) / 2)] # first half
	y_val       = y_notrain[:(len(mols_notrain) / 2)] # first half
	smiles_val  = smiles_notrain[:(len(mols_notrain) / 2)] # first half
	mols_test   = mols_notrain[(len(mols_notrain) / 2):] # second half
	y_test      = y_notrain[(len(mols_notrain) / 2):] # second half
	smiles_test = smiles_notrain[(len(mols_notrain) / 2):] # second half

	# Loop through fingerprint indeces in dense layer
	dense_weights = model.layers[1].get_weights()[0] # weights, not bias
	# for i, dense_weight in enumerate(dense_weights):
	# 	print('at index {}, weight = {}'.format(i, dense_weight))

	# Now report 5 most positive activations:
	sorted_indeces = sorted(range(len(dense_weights)), key = lambda i: dense_weights[i])
	print('Most negatively-activating fingerprint indeces:')
	for j in range(5):
		print('at index {}, weight = {}'.format(sorted_indeces[j], dense_weights[sorted_indeces[j]]))
	print('Most negatively-activating fingerprint indeces:')
	for j in range(1, 6):
		print('at index {}, weight = {}'.format(sorted_indeces[-j], dense_weights[sorted_indeces[-j]]))

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
			test_embeddings_demo(model, config['IO']['data_fpath'], fpath)
			quit(1)
	except KeyError:
		pass

	# Testing activations?
	try:
		if input_to_bool(config['TESTING']['test_activations']):
			# Test current embebddings
			test_activations(model, config['IO']['data_fpath'], fpath)
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
			batch_size = int(config['TRAINING']['batch_size']))
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
