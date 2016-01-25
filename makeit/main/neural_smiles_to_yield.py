from __future__ import print_function
from makeit.utils.parsing import input_to_bool, smiles_to_boolean_fps
from makeit.utils.saving import save_model_history
from keras.models import Graph, Sequential, model_from_json
from keras.layers.core import Dense, Dropout, Activation, Masking
from keras.layers.embeddings import Embedding
from keras.layers.recurrent import LSTM
from keras.preprocessing.sequence import pad_sequences
from keras.optimizers import RMSprop
from keras.utils.visualize_util import plot
import matplotlib.pyplot as plt
import numpy as np
import datetime
import json
import sys
import os

def get_model_fpath():
	'''Returns file path where this model is backed up'''
	fpath_root = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'models')
	fpath = os.path.join(fpath_root, 'neural_smiles_to_yield')
	if flabel:
		return fpath + '_{}'.format(flabel)
	return fpath

def build_model(vocab_size, embedding_size = 100, lstm_size = 32, lr = 0.01):
	'''Generates simple embedding model to use reaction smiles as
	input in order to predict a single-valued output (i.e., yield)'''

	# Base model
	model = Sequential()

	# Add layers
	model.add(Embedding(vocab_size, embedding_size, mask_zero = True))
	print('    model: added Embedding layer ({} -> {})'.format(vocab_size, 
		embedding_size))
	model.add(Dropout(0.2))
	print('    model: added Dropout layer')
	model.add(LSTM(lstm_size))
	print('    model: added LSTM layer ({} -> {})'.format(embedding_size, lstm_size))
	model.add(Dropout(0.2))
	print('    model: added Dropout layer')
	model.add(Dense(1))
	print('    model: added Dense layer ({} -> {})'.format(lstm_size, 1))

	# Compile
	optimizer = RMSprop(lr = lr)
	model.compile(loss = 'mean_squared_error', optimizer = optimizer)

	return model

def save_model(model, data_fpath, hist = None):
	'''Saves NN model object and associated information'''
	fpath = get_model_fpath()

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
	info_fid = open(fpath + '.info', 'w')
	time_now = datetime.datetime.utcnow()
	info_fid.write('{} saved at UTC {}\n\n'.format(fpath, time_now))
	info_fid.write('File details\n------------\n')
	info_fid.write('- data source: {}\n'.format(data_fpath))
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
		data = data[0:10000]

	# Parse data into individual components
	smiles = [x[0] for x in data]
	yields = [x[1] for x in data]
	yields[yields > 1] = 1.0

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
	smiles_train = pad_sequences(smiles_train, maxlen = maxlen)
	smiles_test  = pad_sequences(smiles_test, maxlen = maxlen)

	return (smiles_train, yields_train, smiles_test, yields_test)

def train_model(model, data_fpath, nb_epoch = 3, batch_size = 16):
	'''Trains the model to predict chemical molecular weights'''

	# Get data from helper function
	data = get_data(data_fpath)
	smiles_train = data[0]
	yields_train = data[1]

	# Fit (allows keyboard interrupts in the middle)
	try:
		hist = model.fit(smiles_train, yields_train, 
			nb_epoch = nb_epoch, batch_size = batch_size, verbose = 1)	
	except:
		print('terminated training early due to error')
		hist = None

	return (model, hist)

def test_model(model, data_fpath):
	'''This function evaluates model performance using test data. The output is
	more meaningful than just the loss function value.'''

	# Get data from helper function
	data = get_data(data_fpath)
	smiles_test = data[2]
	yields_test = data[3]

	# Simple evaluation
	score = model.evaluate(smiles_test, yields_test, verbose = 1)
	print('model loss function score: {}'.format(score))

	# Custom evaluation
	yields_predicted = model.predict(smiles_test, verbose = 1)

	# Decode smiles strings
	reverse_tokenizer = {0 : ''}
	for key, val in tokenizer.iteritems():
		reverse_tokenizer[val] = key

	# Save
	fpath = get_model_fpath() + '.test'
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
	plt.title('Parity plot for yield prediction')
	plt.grid(True)
	#plt.axis([0, 1, 0, 1])
	plt.show()

	return score

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print('Usage: {} "tokenizer.json" "data.json" [model_label]'.format(sys.argv[0]))
		print('    tokenizer.json must be a character dictionary')
		print('    data_file.json must be a list of lists, where each' + 
			  ' element is a list of str(reaction_smiles), float(yield))')
		print('    [model_label] is an extra tag to append to the file' + 
			  ' name if a unique identifier is desired')
		quit(1)

	# Load tokenizer
	with open(sys.argv[1], 'rb') as tokenizer_fid:
		tokenizer = json.load(tokenizer_fid)

	# Get model label
	if len(sys.argv) == 4:
		flabel = sys.argv[3]
	else:
		flabel = None

	# See if the model exists in this location already
	fpath = get_model_fpath()
	structure_fpath = fpath + '.json'
	use_old = True
	if os.path.isfile(structure_fpath):
		use_old = raw_input('Use existing model structure [y/n]? ')
		use_old = input_to_bool(use_old)
	else:
		use_old = False
	if use_old:
		# Load model
		with open(structure_fpath, 'r') as structure_fid:
			model = model_from_json(json.load(structure_fid))
			print('...loaded structural information')
	else:
		# Build model
		print('...building model')
		vocab_size = len(tokenizer.keys())
		model = build_model(vocab_size + 1, embedding_size = 100, lstm_size = 100, lr = 0.001)
		print('...built untrained model')

	# See if weights exist in this location already
	weights_fpath = fpath + '.h5'
	if os.path.isfile(weights_fpath):
		use_old = raw_input('Use existing model weights [y/n]? ')
		if input_to_bool(use_old):
			# Use weights
			model.load_weights(weights_fpath)
			print('...loaded weight information')

	# Train model
	print('...training model')
	hist = None
	(model, hist) = train_model(model, sys.argv[2], nb_epoch = 2, batch_size = 500)
	print('...trained model')

	# Save for future
	print('...saving model')
	save_model(model, sys.argv[2], hist = hist)
	print('...saved model')
		
	# Test model
	print('...testing model')
	test_model(model, sys.argv[2])
	print('...tested model')