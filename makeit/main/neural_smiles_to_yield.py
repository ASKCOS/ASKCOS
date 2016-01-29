from __future__ import print_function
from makeit.utils.parsing import input_to_bool, smiles_to_boolean_fps
from makeit.utils.parse_cfg import read_config
from makeit.utils.saving import save_model_history
from keras.models import Graph, Sequential, model_from_json
from keras.layers.core import Dense, Dropout, Activation, Masking
from keras.layers.embeddings import Embedding
from keras.layers.recurrent import LSTM
from keras.preprocessing.sequence import pad_sequences
from keras.optimizers import RMSprop, Adam
from keras.utils.visualize_util import plot
import theano
import matplotlib.pyplot as plt
import numpy as np
import datetime
import json
import sys
import os

def build_model(vocab_size, embedding_size = 100, lstm_size = 100, lr = 0.01, optimizer = 'adam'):
	'''Generates simple embedding model to use reaction smiles as
	input in order to predict a single-valued output (i.e., yield)'''

	# Base model
	model = Sequential()

	# Add layers
	model.add(Embedding(vocab_size, embedding_size, mask_zero = True, init = 'he_uniform'))
	print('    model: added Embedding layer ({} -> {})'.format(vocab_size, 
		embedding_size))
	model.add(LSTM(lstm_size, init = 'he_normal'))
	print('    model: added LSTM layer ({} -> {})'.format(embedding_size, lstm_size))
	model.add(Dropout(0.2))
	print('    model: added Dropout layer')
	model.add(Dense(1, init = 'zero'))
	print('    model: added Dense layer ({} -> {})'.format(lstm_size, 1))

	# Compile
	if optimizer == 'adam':
		optimizer = Adam(lr = lr)
	elif optimizer == 'rmsprop':
		optimizer = RMSProp(lr = lr)
	else:
		print('Can only handle adam or rmsprop optimizers currently')
		quit(1)

	model.compile(loss = 'mean_squared_error', optimizer = optimizer)

	return model

def save_model(model, fpath = '', config = {}, hist = None, tstamp = ''):
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
	if hist:
		save_model_history(hist, fpath + '.hist')
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
	yields = np.array([x[1] for x in data], dtype = np.float32)
	yields[yields > 1.0] = 1.0

	# Vectorize chemical names according to tokenizer
	smiles = [[tokenizer[x] for x in one_smiles] for one_smiles in smiles]
	maxlen = max([len(x) for x in smiles])

	# Create training/development split
	division = int(len(data) * training_ratio)
	smiles_train = smiles[:division]
	smiles_test  = smiles[division:]
	yields_train = yields[:division]
	yields_test  = yields[division:]

	# Pad to get uniform lengths (req. for training in same batch)
	smiles_train = pad_sequences(smiles_train, maxlen = maxlen).astype(np.uint8)
	smiles_test  = pad_sequences(smiles_test, maxlen = maxlen).astype(np.uint8)

	return (smiles_train, yields_train, smiles_test, yields_test)

def train_model(model, data_fpath = '', nb_epoch = 0, batch_size = 1):
	'''Trains the model to predict chemical molecular weights'''

	# Get data from helper function
	data = get_data(data_fpath)
	smiles_train = data[0]
	yields_train = data[1]

	# Fit (allows keyboard interrupts in the middle)
	try:
		hist = model.fit(smiles_train, yields_train, 
			nb_epoch = nb_epoch, 
			batch_size = batch_size, 
			verbose = 1)	
	except KeyboardInterrupt:
		print('terminated training early (intentionally)')
		hist = None

	return (model, hist)

def test_model(model, data_fpath, fpath, tstamp = '', batch_size = 128):
	'''This function evaluates model performance using test data. The output is
	more meaningful than just the loss function value.'''
	test_fpath = fpath + ' ' + tstamp

	# Get data from helper function
	data = get_data(data_fpath)
	smiles_train = data[0]
	yields_train = data[1]
	smiles_test = data[2]
	yields_test = data[3]

	# Simple evaluation
	score = model.evaluate(smiles_test, yields_test, 
		batch_size = batch_size,
		verbose = 1)
	print('model loss function score: {}'.format(score))

	# Custom evaluation
	yields_predicted = model.predict(smiles_test, 
		batch_size = batch_size,
		verbose = 1)
	yields_train_predicted = model.predict(smiles_train, 
		batch_size = batch_size,
		verbose = 1)

	# Decode smiles strings
	reverse_tokenizer = {0 : ''}
	for key, val in tokenizer.iteritems():
		reverse_tokenizer[val] = key

	# Save
	with open(fpath, 'w') as fid:
		time_now = datetime.datetime.utcnow()
		fid.write('-- tested at UTC {}\n\n'.format(fpath, time_now))		
		fid.write('test entry\treaction\tactual yield\tpredicted yield\n')
		for i in range(len(smiles_test)):
			fid.write('{}\t{}\t{}\t{}\n'.format(i, 
				''.join([reverse_tokenizer[x] for x in smiles_test[i]]), 
				yields_test[i], yields_predicted[i, 0]))

	# Create parity plot
	plt.scatter(yields_test, yields_predicted, alpha = 0.2)
	plt.xlabel('Actual yield')
	plt.ylabel('Predicted yield')
	plt.title('Parity plot for yield prediction (test set)')
	plt.grid(True)
	plt.axis([0, 1, 0, 1])
	plt.savefig(test_fpath + ' test.png', bbox_inches = 'tight')
	plt.clf()

	# Look at training data onw
	plt.scatter(yields_train, yields_train_predicted, alpha = 0.2)
	plt.xlabel('Actual yield')
	plt.ylabel('Predicted yield')
	plt.title('Parity plot for yield prediction (train set)')
	plt.grid(True)
	plt.axis([0, 1, 0, 1])
	plt.savefig(test_fpath + ' train.png', bbox_inches = 'tight')
	plt.clf()

	return score

def test_embeddings_demo(model, data_fpath, ref_index = 0):
	'''This function tests dense embeddings of reactions and tries to find the
	most similar one to a test example'''

	# Get training data from helper function
	data = get_data(data_fpath)
	smiles_train = data[0]
	smiles_test = data[2]

	# Decode smiles strings
	reverse_tokenizer = {0 : ''}
	for key, val in tokenizer.iteritems():
		reverse_tokenizer[val] = key

	# Define function to test embedding
	tf = theano.function([model.layers[0].input], 
		model.layers[1].get_output(train = False))

	# Look for reaction most similar to test reaction one:
	ref_embedding = tf([smiles_test[ref_index]])[0]
	print('shape of model.layers[1] output for comparison: {}'.format(ref_embedding.shape))

	# For debugging
	np.set_printoptions(threshold = 'nan')
	print(ref_embedding)
	
	# Define vector-angle function
	def vector_angle(embedding):
		angle = np.arccos(np.dot(embedding, ref_embedding) / 
			np.linalg.norm(embedding) / 
			np.linalg.norm(ref_embedding))
		return angle

	current_closest = [[0]]
	current_minangle  = 99
	nb = 2000
	i = 0
	batches = [smiles_train[x:(x + nb)] for x in range(0, len(smiles_train), nb)]
	N = len(batches)
	for batch in batches:
		embeddings = tf(batch)
		i = i + 1
		print('{}/{}'.format(i, N))
		scores = np.apply_along_axis(vector_angle, axis = 1, arr = embeddings)
		minangle_index = np.argmin(scores)
		minangle = scores[minangle_index]
		if minangle < current_minangle:
			current_closest = batch[minangle_index]
			current_minangle = minangle
			print('...found a better match, angle = {}'.format(current_minangle))
			# For debugging
			np.set_printoptions(threshold='nan')
			print(embeddings[minangle_index])
	print('---results---')
	print('Reference reaction: {}'.format(''.join([reverse_tokenizer[x] for x in smiles_test[ref_index]])))
	print('Most similar from training: {}'.format(''.join([reverse_tokenizer[x] for x in current_closest])))
	print(' (angle = {})'.format(current_minangle))

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

	# Load tokenizer
	with open(config['IO']['tokenizer_fpath'], 'rb') as tokenizer_fid:
		tokenizer = json.load(tokenizer_fid)

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
			model = build_model(len(tokenizer.keys()) + 1, 
				embedding_size = int(config['ARCHITECTURE']['embedding_size']), 
				lstm_size = int(config['ARCHITECTURE']['lstm_size']), 
				lr = float(config['ARCHITECTURE']['lr']),
				optimizer = config['ARCHITECTURE']['optimizer'])
			print('...built untrained model')
		except KeyError:
			print('Must specify embedding_size, lstm_size, and lr in ARCHITECTURE in config')
			quit(1)
		except ValueError:
			print('embedding_size, lstm_size, and/or lr must be suitable numbers')
			quit(1)

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
			test_embeddings_demo(model, config['IO']['data_fpath'],
				ref_index = int(config['TESTING']['test_index']))
			quit(1)
	except KeyError:
		pass

	# Train model
	hist = None
	try:
		print('...training model')
		model.optimizer.lr.set_value(float(config['TRAINING']['lr']))
		(model, hist) = train_model(model, 
			data_fpath = config['IO']['data_fpath'], 
			nb_epoch = int(config['TRAINING']['nb_epoch']), 
			batch_size = int(config['TRAINING']['batch_size']))
		print('...trained model')
	except KeyError:
		print('Must specify data_fpath in IO and nb_epoch and batch_size in TRAINING in config')
		quit(1)
	except ValueError:
		print('nb_epoch and batch_size must be integers')
		quit(1)
	except KeyboardInterrupt:
		pass

	# Get the current time
	tstamp = datetime.datetime.utcnow().strftime('%m-%d-%Y_%H-%M')

	# Save for future
	print('...saving model')
	save_model(model, 
		fpath = fpath,
		config = config, 
		hist = hist, 
		tstamp = tstamp)
	print('...saved model')
		
	# Test model
	print('...testing model')
	test_model(model, config['IO']['data_fpath'],
		fpath,
		tstamp = tstamp,
		batch_size = int(config['TRAINING']['batch_size']))
	print('...tested model')
