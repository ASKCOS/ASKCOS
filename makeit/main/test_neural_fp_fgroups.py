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
import theano # for looking at computational graph
import rdkit.Chem as Chem
import matplotlib.pyplot as plt
import numpy as np
import datetime
import json
import sys
import os

def build_model(embedding_size = 100, lr = 0.01, optimizer = 'adam'):
	'''Generates simple embedding model to use reaction smiles as
	input in order to predict a single-valued output (i.e., yield)'''

	# Base model
	model = Sequential()

	# Add layers
	model.add(GraphFP(embedding_size, sizeAttributeVector() - 1, depth = 2))
	print('    model: added GraphFP layer ({} -> {})'.format('mol', embedding_size))
	model.add(Activation('tanh'))
	print('    model: added tanh layer ({} -> {})'.format(embedding_size, embedding_size))

	# Compile
	if optimizer == 'adam':
		optimizer = Adam(lr = lr)
	elif optimizer == 'rmsprop':
		optimizer = RMSProp(lr = lr)
	else:
		print('Can only handle adam or rmsprop optimizers currently')
		quit(1)

	model.compile(loss = 'mse', optimizer = optimizer)

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
	print('...loading data')
	with open(data_fpath, 'r') as data_fid:
		data = json.load(data_fid)
		
	# Truncate if necessary
	try:
		truncate_to = int(config['TRAINING']['truncate_to'])
		data = data[0:truncate_to]
	except:
		pass

	# Get new training_ratio if possible
	try:
		temp = float(config['TRAINING']['ratio'])
		training_ratio = temp
	except:
		pass

	# Parse data into individual components
	smiles = [x[0] for x in data]
	mols = [molToGraph(Chem.MolFromSmiles(smile)).dump_as_tensor() for smile in smiles]
	fgroups = [np.array(x[1]) > 0 for x in data] # change count to flag

	# Create training/development split
	division = int(len(data) * training_ratio)
	mols_train = mols[:division]
	mols_test  = mols[division:]
	fgroups_train = fgroups[:division]
	fgroups_test  = fgroups[division:]
	smiles_train = smiles[:division]
	smiles_test = smiles[division:]

	return (mols_train, fgroups_train, mols_test, fgroups_test, training_ratio, smiles_train, smiles_test)

def train_model(model, data_fpath = '', nb_epoch = 0, batch_size = 1):
	'''Trains the model'''

	# Get data from helper function
	data = get_data(data_fpath)
	mols_train = data[0]
	fgroups_train = data[1]
	mols_test = data[2]
	fgroups_test = data[3]
	training_ratio = data[4]
	# Note: will rejoin data and use validation_split

	# Fit (allows keyboard interrupts in the middle)
	# Because molecular graph tensors are different sizes based on N_atoms, can only do one at a time
	# (alternative is to pad with zeros and try to add some masking feature to GraphFP)
	try:
		loss = []
		val_loss = []
		for i in range(nb_epoch):
			print('Epoch {}/{}'.format(i + 1, nb_epoch))
			this_loss = []
			this_val_loss = []
			
			print('Training...')
			for j in range(len(mols_train)):
				single_mol_as_array = np.array(mols_train[j:j+1])
				single_fp_as_array = np.array(fgroups_train[j:j+1])
				sloss = model.train_on_batch(single_mol_as_array, single_fp_as_array, accuracy = False)
				# print('single loss: {}'.format(sloss))
				this_loss.append(sloss)
			
			print('Testing...')
			for j in range(len(mols_test)):
				single_mol_as_array = np.array(mols_test[j:j+1])
				spred = model.predict_on_batch(single_mol_as_array)
				sloss = np.mean((spred[0] - fgroups_test[j]) ** 2)
				# print('single loss: {}'.format(sloss))
				this_val_loss.append(sloss)
			
			loss.append(np.mean(this_loss))
			val_loss.append(np.mean(this_val_loss))
			print('loss: {}\tval_loss: {}'.format(loss[i], val_loss[i]))

		# Write history

	except KeyboardInterrupt:
		print('terminated training early (intentionally)')

	return (model, loss, val_loss)

def test_model(model, data_fpath, fpath, tstamp = '', batch_size = 128):
	'''This function evaluates model performance using test data. The output is
	more meaningful than just the loss function value.'''
	# test_fpath = fpath + ' ' + tstamp

	# # Get data from helper function
	# data = get_data(data_fpath)
	# smiles_train = data[0]
	# yields_train = data[1]
	# smiles_test = data[2]
	# yields_test = data[3]

	# # Custom evaluation
	# yields_predicted = model.predict(smiles_test, 
	# 	batch_size = batch_size,
	# 	verbose = 1)
	# yields_train_predicted = model.predict(smiles_train, 
	# 	batch_size = batch_size,
	# 	verbose = 1)

	# # Save
	# with open(fpath, 'w') as fid:
	# 	time_now = datetime.datetime.utcnow()
	# 	fid.write('-- tested at UTC {}\n\n'.format(fpath, time_now))		
	# 	fid.write('test entry\treaction\tactual yield\tpredicted yield\n')
	# 	for i in range(len(smiles_test)):
	# 		fid.write('{}\t{}\t{}\t{}\n'.format(i, 
	# 			''.join([reverse_tokenizer[x] for x in smiles_test[i]]), 
	# 			yields_test[i], yields_predicted[i, 0]))

	# # Create parity plot
	# plt.scatter(yields_test, yields_predicted, alpha = 0.2)
	# plt.xlabel('Actual yield')
	# plt.ylabel('Predicted yield')
	# plt.title('Parity plot for yield prediction (test set)')
	# plt.grid(True)
	# plt.axis([0, 1, 0, 1])
	# plt.savefig(test_fpath + ' test.png', bbox_inches = 'tight')
	# plt.clf()

	# # Look at training data onw
	# plt.scatter(yields_train, yields_train_predicted, alpha = 0.2)
	# plt.xlabel('Actual yield')
	# plt.ylabel('Predicted yield')
	# plt.title('Parity plot for yield prediction (train set)')
	# plt.grid(True)
	# plt.axis([0, 1, 0, 1])
	# plt.savefig(test_fpath + ' train.png', bbox_inches = 'tight')
	# plt.clf()

	return

def test_embeddings_demo(model, data_fpath):
	'''This function tests dense embeddings of reactions and tries to find the
	most similar one to a test example'''

	# Get data from helper function
	data = get_data(data_fpath)
	mols_train = data[0]
	fgroups_train = data[1]
	mols_test = data[2]
	fgroups_test = data[3]
	training_ratio = data[4]
	smiles_train = data[5]
	smiles_test = data[6]

	print('-------------')
	print('DEMONSTRATION')
	print('-------------')
	for j in range(len(mols_train)):
		print('original smiles: {}'.format(smiles_train[j]))
		embedding = model.predict_on_batch(np.array(mols_train[j:j+1]))[0][0]
		print('{}\t{}\t{}'.format('idx', 'embed', 'target'))
		for k in range(len(embedding)):
			print('{}\t{}\t{}'.format(k, round(embedding[k], 2), int(fgroups_train[j][k])))

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
			model = build_model(embedding_size = int(config['ARCHITECTURE']['embedding_size']))
			print('...built untrained model')
		except KeyError:
			print('Must specify embedding_size in ARCHITECTURE in config')
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
