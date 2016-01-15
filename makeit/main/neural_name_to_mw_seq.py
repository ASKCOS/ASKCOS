from __future__ import print_function
from makeit.utils.parsing import SplitChemicalName, input_to_bool, sequences_to_texts
from makeit.utils.saving import save_model_history
from keras.models import Sequential, model_from_json
from keras.layers.core import Dense, Dropout, Activation
from keras.layers.embeddings import Embedding
from keras.preprocessing.sequence import pad_sequences
from keras.layers.recurrent import LSTM
from keras.optimizers import SGD
from keras.utils.visualize_util import plot
import datetime
import cPickle
import json
import sys
import os

def get_model_fpath():
	'''Returns file path where this model is backed up'''
	fpath_root = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'models')
	fpath = os.path.join(fpath_root, 'neural_name_to_mw_seq')
	if flabel:
		return fpath + '_{}'.format(flabel)
	return fpath

def build_model(vocab_size, embedding_size = 100, lstm_size = 16):
	'''Generates simple embedding model to use tokenized chemical name as
	input in order to predict a single-valued output (i.e., MW)'''

	# Base model
	model = Sequential()

	# Add layers
	model.add(Embedding(vocab_size, embedding_size, mask_zero = True))
	print('    model: added Embedding layer ({} -> {})'.format(vocab_size, 
		embedding_size))
	# model.add(Dropout(0.5))
	# print('    model: added Dropout layer')
	model.add(LSTM(lstm_size))
	print('    model: added LSTM layer ({} -> {})'.format(embedding_size, lstm_size))
	# model.add(Activation('linear'))
	# print('    model: added Activation(linear) layer')
	model.add(Dropout(0.5))
	print('    model: added Dropout layer')
	model.add(Dense(1))
	print('    model: added Dense layer ({} -> {})'.format(lstm_size, 1))

	# Compile
	model.compile(loss = 'mean_squared_error', optimizer = 'rmsprop')

	return model

def save_model(model, tokenizer_fpath, data_fpath, hist = None):
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
	info_fid.write('- tokenizer source: {}\n'.format(tokenizer_fpath))
	info_fid.write('- data source: {}\n'.format(data_fpath))
	info_fid.close()

	print('...saved model to {}.[json, h5, png, info]'.format(fpath))
	return True

def get_data(data_fpath, training_ratio = 0.9, maxlen = 30):
	'''This is a helper script to read the data .json object and return
	the training and test data sets separately. This is to allow for an
	already-trained model to be evaluated using the test data (i.e., which
	we know it hasn't seen before)'''

	# Load data from json file
	with open(data_fpath, 'r') as data_fid:
		data = json.load(data_fid)

	# Parse data into individual components
	names = [' '.join(SplitChemicalName(str(x[0]))) for x in data] # from unicode
	mws   = [x[1] for x in data]

	# Vectorize chemical names according to tokenizer
	names = tokenizer.texts_to_sequences(names)

	# Create training/development split
	division = int(len(data) * training_ratio)
	names_train = names[:division]
	names_test  = names[division:]
	mws_train = mws[:division]
	mws_test  = mws[division:]

	# Pad to get uniform lengths (req. for training in same batch)
	names_train = pad_sequences(names_train, maxlen = maxlen)
	names_test  = pad_sequences(names_test, maxlen = maxlen)

	return (names_train, mws_train, names_test, mws_test)

def train_model(model, tokenizer, data_fpath, nb_epoch = 3, batch_size = 16):
	'''Trains the model to predict chemical molecular weights'''

	# Get data from helper function
	data = get_data(data_fpath)
	names_train = data[0]
	mws_train = data[1]
	# Fit
	hist = model.fit(names_train, mws_train, nb_epoch = nb_epoch, 
			batch_size = batch_size, verbose = 1)	

	return (model, hist)

def test_model(model, tokenizer, data_fpath):
	'''This function evaluates model performance using test data. The output is
	more meaningful than just the loss function value.'''

	# Get data from helper function
	data = get_data(data_fpath)
	names_test = data[2]
	mws_test = data[3]

	# Simple evaluation
	score = model.evaluate(names_test, mws_test, verbose = 1)
	print('model loss function score: {}'.format(score))

	# Custom evaluation
	mws_predicted = model.predict(names_test, verbose = 1)

	names_test_decoded = sequences_to_texts(tokenizer, names_test)
	for i in range(5):
		print('Test entry {}'.format(i))
		print('  decoded name: {}'.format(names_test_decoded[i]))
		print('  actual mw:    {}'.format(mws_test[i]))
		print('  predicted mw: {}'.format(mws_predicted[i, 0]))

	# Save
	fpath = get_model_fpath() + '.test'
	with open(fpath, 'w') as fid:
		fid.write('Detailed model test:\n')
		fid.write('test entry\tdecoded name\tactual mw\tpredicted mw\n')
		for i in range(len(names_test)):
			fid.write('{}\t{}\t{}\t{}\n'.format(i, names_test_decoded[i], 
					mws_test[i], mws_predicted[i, 0]))

	return score

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print('Usage: {} "tokenizer.cpickle" "data.json" [model_label]'.format(sys.argv[0]))
		print('    tokenizer.cpickle must be an already-fit tokenizer')
		print('    data_file.json must be a list of lists, where each' + 
			  ' element is a list of str(name) and float(mw)')
		print('    [model_label] is an extra tag to append to the file' + 
			  ' name if a unique identifier is desired')
		quit(1)

	# Get model label
	if len(sys.argv) == 4:
		flabel = sys.argv[3]
	else:
		flabel = None

	# Load tokenizer
	with open(sys.argv[1], 'rb') as tokenizer_fid:
		tokenizer = cPickle.load(tokenizer_fid)

	# See if the model exists in this location already
	fpath = get_model_fpath()
	structure_fpath = fpath + '.json'
	use_old = True
	if os.path.isfile(structure_fpath):
		use_old = raw_input('Use existing model structure [y/n]? ')
		use_old = input_to_bool(use_old)
	if use_old:
		# Load model
		with open(structure_fpath, 'r') as structure_fid:
			model = model_from_json(json.load(structure_fid))
			print('...loaded structural information')
	else:
		# Build model
		print('...building model')
		if tokenizer.nb_words:
			vocab_size = min([tokenizer.nb_words, len(tokenizer.word_counts)])
		else:
			vocab_size = len(tokenizer.word_counts)
		model = build_model(vocab_size)
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
	(model, hist) = train_model(model, tokenizer, sys.argv[2], nb_epoch = 1, batch_size = 64)
	print('...trained model')

	# Save for future
	print('...saving model')
	save_model(model, sys.argv[1], sys.argv[2], hist = hist)
	print('...saved model')
		
	# Test model
	print('...testing model')
	test_model(model, tokenizer, sys.argv[2])
	print('...tested model')