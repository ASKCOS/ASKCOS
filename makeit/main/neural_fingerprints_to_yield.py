from __future__ import print_function
from makeit.utils.parsing import input_to_bool, smiles_to_boolean_fps
from makeit.utils.saving import save_model_history
from keras.models import Graph, model_from_json
from keras.layers.core import Dense, Dropout, Activation
from keras.layers.embeddings import Embedding
from keras.layers.recurrent import LSTM
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
	fpath = os.path.join(fpath_root, 'neural_fingerprints_to_yield')
	if flabel:
		return fpath + '_{}'.format(flabel)
	return fpath

def build_model(fp_size = 2048, embedding_size = 100, dropout = 0.75):
	'''Uses fingerprints for two chemicals to predict the yield of a reaction
	between them, ignoring the chemical species of that product'''

	# Base model
	model = Graph()

	# Inputs
	model.add_input(name = 'fpA', input_shape = (fp_size, ), dtype = 'bool')
	model.add_input(name = 'fpB', input_shape = (fp_size, ), dtype = 'bool')

	# Embebdding
	embedding = Embedding(output_dim = embedding_size, input_dim = fp_size, 
		init = 'glorot_normal')
	model.add_shared_node(embedding, 'embed', inputs = ['fpA', 'fpB'], 
		outputs = ['emA', 'emB'])
	model.add_shared_node(Dropout(dropout), 'dropout', inputs = ['emA', 'emB'], 
		outputs = ['drA', 'drB'])
	model.add_shared_node(LSTM(output_dim = embedding_size, activation = 'tanh', 
		return_sequences = True), 'lstm', inputs = ['drA', 'drB'], 
		outputs = ['lstmA', 'lstmB'])

	# Merge
	model.add_node(Dropout(dropout), name = 'drM', inputs = ['lstmA', 'lstmB'], 
		merge_mode = 'concat')
	model.add_node(LSTM(output_dim = 2 * embedding_size, activation = 'tanh', 
		return_sequences = False), 'lstmM', input = 'drM')
	model.add_node(Dropout(dropout), name = 'drmM2', input = 'lstmM')
	model.add_node(Dense(1, activation = 'linear'), name = 'denseM',
		input = 'drmM2')
	model.add_output(name = 'yield', input = 'denseM')

	# Compile
	print('...compiling model')
	optimizer = RMSprop(lr = 0.1)
	model.compile(loss = {'yield' : 'mse'}, optimizer = optimizer)

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
		data = data[0:1000]

	# Parse data into individual components
	print('...calculating fingerprints')
	fpA = smiles_to_boolean_fps([x[0] for x in data])
	fpB = smiles_to_boolean_fps([x[1] for x in data])
	yields   = np.asarray([x[2] for x in data])

	# Create training/development split
	division = int(len(data) * training_ratio)
	fpA_train = fpA[:division]
	fpA_test  = fpA[division:]
	fpB_train = fpB[:division]
	fpB_test  = fpB[division:]
	yields_train = yields[:division]
	yields_test  = yields[division:]
	print('...converted data')

	return (fpA_train, fpB_train, yields_train, 
		    fpA_test, fpB_test, yields_test)

def train_model(model, data_fpath, nb_epoch = 3, batch_size = 16):
	'''Trains the model to predict chemical molecular weights'''

	# Get data from helper function
	data = get_data(data_fpath)
	fpA_train = data[0]
	fpB_train = data[1]
	yields_train = data[2]

	# Fit
	hist = model.fit({'fpA' : fpA_train, 'fpB' : fpB_train, 'yield' : yields_train}, 
		nb_epoch = nb_epoch, batch_size = batch_size, verbose = 1)	

	return (model, hist)

def test_model(model, data_fpath):
	'''This function evaluates model performance using test data. The output is
	more meaningful than just the loss function value.'''

	# Get data from helper function
	data = get_data(data_fpath)
	fpA_test = data[3]
	fpB_test = data[4]
	yields_test = data[5]

	# Simple evaluation
	score = model.evaluate({'fpA' : fpA_test, 'fpB' : fpB_test, 'yield' : yields_test}, 
		verbose = 1)
	print('model loss function score: {}'.format(score))

	# Custom evaluation
	yields_predicted = model.predict({'fpA' : fpA_test, 'fpB' : fpB_test}, 
		verbose = 1)['yield']

	# Save
	fpath = get_model_fpath() + '.test'
	with open(fpath, 'w') as fid:
		fid.write('Detailed model test:\n')
		fid.write('test entry\tactual yield\tpredicted yield\n')
		for i in range(len(fpA_test)):
			fid.write('{}\t{}\t{}\t{}\n'.format(i, yields_test[i], 
					yields_predicted[i]))

	# Create parity plot
	plt.scatter(yields_test, yields_predicted, alpha = 0.2)
	plt.xlabel('Actual yield')
	plt.ylabel('Predicted yield')
	plt.title('Parity plot for yield prediction')
	plt.grid(True)
	plt.axis([0, 1, 0, 1])
	plt.show()

	return score

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print('Usage: {} "data.json" [model_label]'.format(sys.argv[0]))
		print('    data_file.json must be a list of lists, where each' + 
			  ' element is a list of str(A_smiles), str(B_smiles) and float(yield))')
		print('    [model_label] is an extra tag to append to the file' + 
			  ' name if a unique identifier is desired')
		quit(1)

	# get_data(sys.argv[1])
	# quit(1)

	# Get model label
	if len(sys.argv) == 3:
		flabel = sys.argv[2]
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
		model = build_model(fp_size = 2048, embedding_size = 100, dropout = 0.75)
		print('...built untrained model')
		# Save for future
		print('...saving model')
		save_model(model, sys.argv[1])
		print('...saved model')

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
	(model, hist) = train_model(model, sys.argv[1], nb_epoch = 1, batch_size = 16)
	print('...trained model')

	# Save for future
	print('...saving model')
	save_model(model, sys.argv[1], hist = hist)
	print('...saved model')
		
	# Test model
	print('...testing model')
	test_model(model, sys.argv[1])
	print('...tested model')