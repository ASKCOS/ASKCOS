from __future__ import print_function
from makeit.utils.parsing import SplitChemicalName 
from keras.models import Sequential
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
	return os.path.join(fpath_root, 'neural_name_to_mw')

def build_model(vocab_size, embedding_size = 100):
	'''Generates simple embedding model to use tokenized chemical name as
	input in order to predict a single-valued output (i.e., MW)'''

	# Base model
	model = Sequential()

	# Add layers
	model.add(Embedding(vocab_size, embedding_size, 
		init = 'uniform', mask_zero = True)
	print('    model: added Embedding layer ({} -> {})'.format(vocab_size, 
		embedding_size))
	model.add(Dropout(0.5))
	print('    model: added Dropout layer')
	model.add(LSTM(output_dim = 1))
	print('    model: added LSTM layer ({} -> {})'.format(embedding_size, 1))

	# Define optimizer
	sgd = SGD(lr = 0.1, decay = 1e-6, momentum = 0.9)

	# Compile
	model.compile(loss = 'mean_squared_error', optimizer = sgd)

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
		hist_fid = open(fpath + '.hist', 'w')
		print(hist.history, file = hist_fid)

	# Write to info file
	info_fid = open(fpath + '.info', 'w')
	time_now = datetime.datetime.utcnow()
	info_fid.write('{} generated at UTC {}\n\n'.format(fpath, time_now))
	info_fid.write('File details\n------------\n')
	info_fid.write('- tokenizer source: {}\n'.format(tokenizer_fpath))
	info_fid.write('- data source: {}\n'.format(data_fpath))
	info_fid.close()

	print('...saved model to {}.[json, h5, png, info]'.format(fpath))
	return True

def train_model(model, tokenizer, data_fpath):
	'''Trains the model to predict chemical molecular weights'''

	# Load data from json file
	with open(data_fpath, 'r') as data_fid:
		data = json.load(data_fid)

	# Parse data into individual components
	names = [str(x[0]) for x in data] # from unicode
	mws   = [x[1] for x in data]

	# Create training/development split
	division = int(len(data) * 0.9)
	names_train = names[:division]
	names_test  = names[division:]
	mws_train = mws[:division]
	mws_test  = mws[division:]

	# Vectorize chemical names according to tokenizer
	names_train = tokenizer.texts_to_sequences(names_train)
	names_test  = tokenizer.texts_to_sequences(names_test)
	print('names_train shape: {}'.format(names_train.shape))
	print('names_test  shape: {}'.format(names_test.shape))

	# Pad to get uniform lengths (req. for training in same batch)
	names_train = pad_sequences(names_train, maxlen = 30)
	names_test  = pad_sequences(names_test, maxlen = 30)
	print('names_train shape: {}'.format(names_train.shape))
	print('names_test  shape: {}'.format(names_test.shape))

	#return (model, -1)

	hist = model.fit(names_train, mws_train, nb_epoch = 5, batch_size = 4, 
	 	      verbose = 1, show_accuracy = True)
	score = model.evaluate(names_test, mws_test, batch_size = 4, 
			  verbose = 1, show_accuracy = True)

	return (model, hist, score)


if __name__ == '__main__':
	if len(sys.argv) < 3:
		print('Usage: {} "tokenizer.cpickle" "data.json"'.format(sys.argv[0]))
		print('    tokenizer.cpickle must be an already-fit tokenizer')
		print('    data_file.json must be a list of lists, where each' + 
			  ' element is a list of str(name) and float(mw)')
		quit(1)

	# Load tokenizer
	with open(sys.argv[1], 'rb') as tokenizer_fid:
		tokenizer = cPickle.load(tokenizer_fid)

	# Build model
	print('...building model')
	if tokenizer.nb_words:
		vocab_size = min([tokenizer.nb_words, len(tokenizer.word_counts)])
	else:
		vocab_size = len(tokenizer.word_counts)
	model = build_model(vocab_size)
	print('...built untrained model')

	# Train model
	print('...training model')
	(model, hist, score) = train_model(model, tokenizer, sys.argv[2])
	print('...trained model, evaluated score = {}'.format(score))

	# Save for future
	print('...saving model')
	save_model(model, sys.argv[1], sys.argv[2], hist = hist)
	print('...saved model')
		

