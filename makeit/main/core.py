from __future__ import print_function
from makeit.utils.saving import save_model_history, save_model_history_manual
from makeit.utils.GraphEmbedding import GraphFP
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Masking
from keras.callbacks import LearningRateScheduler, EarlyStopping
from keras.regularizers import l2
from keras.optimizers import RMSprop, Adam
import numpy as np
import datetime
import json
import sys
import os


def build_model(embedding_size = 100, lr = 0.01, optimizer = 'adam', depth = 2, 
	scale_output = 0.05, padding = True, inner_reg_l2 = 0.0, output_reg_l2 = 0.0,
	hidden = 0, loss = 'mse'):
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
		model.add(Dense(hidden, activation = 'tanh'))
		print('    model: added Dense layer (-> {})'.format(hidden))
	model.add(Dense(1))
	print('    model: added Dense layer (-> {})'.format(1))

	# Compile
	if optimizer == 'adam':
		optimizer = Adam(lr = lr)
	elif optimizer == 'rmsprop':
		optimizer = RMSProp(lr = lr)
	else:
		print('Can only handle adam or rmsprop optimizers currently')
		quit(1)

	print('compiling...',)
	model.compile(loss = loss, optimizer = optimizer)
	print('done')

	return model

def save_model(model, loss, val_loss, fpath = '', config = {}, tstamp = ''):
	'''Saves NN model object and associated information.

	inputs:
		model - a Keras model
		loss - list of training losses 
		val_loss - list of validation losses
		fpath - root filepath to save everything to (with .json, h5, png, info 
		config - the configuration dictionary that defined this model 
		tstamp - current timestamp to log in info file'''

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


def train_model(model, get_data, nb_epoch = 0, batch_size = 1, lr_func = '0.01', patience = 10):
	'''Trains the model.

	inputs:
		model - a Keras model
		get_data - the function which returns three dictionaries for training,
				validation, and testing separately
		nb_epoch - number of epochs to train for
		batch_size - batch_size to use on the data. This must agree with what was
				specified for get_data (i.e., if tensors are padded or not)
		lr_func - string which is evaluated with 'epoch' to produce the learning 
				rate at each epoch 
		patience - number of epochs to wait when no progress is being made in 
				the validation loss

	outputs:
		model - a trained Keras model
		loss - list of training losses corresponding to each epoch 
		val_loss - list of validation losses corresponding to each epoch'''

	# Get data from helper function
	(train, val, test) = get_data()
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
				print('mse loss: {}\tmse val_loss: {}'.format(loss[i], val_loss[i]))

				# Check progress
				if np.mean(this_val_loss) < prev_best_val_loss:
					wait = 0
					prev_best_val_loss = np.mean(this_val_loss)
				else:
					wait = wait + 1
					print('{} epochs without val_loss progress'.format(wait))
					if wait == patience:
						print('stopping early!')
						break

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
