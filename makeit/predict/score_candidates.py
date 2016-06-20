# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import numpy as np
import os
import sys
from keras import backend as K 
from keras.models import Sequential
from keras.layers import Dense, Activation
from keras.layers.core import Flatten, Permute, Reshape, Dropout
from keras.layers.wrappers import TimeDistributed
from keras.optimizers import *
from keras.layers.convolutional import Convolution1D, Convolution2D
from keras.regularizers import l2
from makeit.predict.preprocess_candidates import *
import theano.tensor as T
from scipy.sparse import coo_matrix
import cPickle as pickle

def build(N_e = 2048, N_h1 = 100, N_h2 = 50, N_c = 10000, inner_act = 'tanh',
		dr = 0.5):
	'''
	Builds the feed forward model.

	N_e:  size of the fixed reaction embedding
	N_h1: number of hidden nodes in first layer
	N_h2: number of hidden nodes in second layer
	N_c:  maximum number of candidates for each sample (padded length)
	inner_act: activation function 
	'''

	model = Sequential()

	# Time distributed dense
	model.add(TimeDistributed(
		Dense(N_h1, activation = inner_act, W_regularizer = l2(0.001)), 
		input_shape = (N_c, N_e))
	)
	model.add(Dropout(dr))

	# Time distributed dense
	if N_h2 > 0:
		model.add(TimeDistributed(Dense(N_h2, activation = inner_act)))
		model.add(Dropout(dr))

	# Time distributed dense
	model.add(TimeDistributed(Dense(1, activation = 'linear')))

	# Flatten to 1 x N_c
	model.add(Flatten())

	# Softmax
	model.add(Activation('softmax'))

	# Now compile
	sgd = SGD(lr = 0.001, decay = 1e-4, momentum = 0.9)
	# model.compile(loss = 'binary_crossentropy', optimizer = sgd, 
	# 	metrics = ['binary_accuracy'])
	model.compile(loss = 'categorical_crossentropy', optimizer = sgd, 
		metrics = ['accuracy'])

	return model

def train(model, x, y):

	nb_epoch = 100
	batch_size = 40

	try:
		hist = model.fit(x, y, 
			nb_epoch = nb_epoch, 
			batch_size = batch_size, 
			validation_split = 0.2,
			verbose = 1)
	except KeyboardInterrupt:
		print('Stopped training early!')
		return None

	return hist


if __name__ == '__main__':

	if len(sys.argv) >= 2:
		tag = '_' + sys.argv[1].strip()
	else:
		tag = ''

	# Data must be pre-padded
	with open('x_coo_singleonly.dat', 'rb') as infile:
		x = pickle.load(infile)
		x = [np.array(coo_matrix((z[0], (z[1], z[2])), shape = z[3]).todense()) for z in x]
		x = np.transpose(np.array(x), (0, 2, 1))
	with open('y_coo_singleonly.dat', 'rb') as infile:
		y = pickle.load(infile)
		y = np.array(coo_matrix((y[0], (y[1], y[2])), shape = y[3]).todense())

	# Build and train
	N_c = len(y[0,:])
	print('Padded number of candidates for each (N_c): {}'.format(N_c))
	print('Mean number of true predictions for each:   {}'.format(np.mean(np.sum(y, axis = 1))))

	model = build(N_c = N_c)
	hist = train(model, x, y)
	model.save_weights('weights{}.h5'.format(tag))

	# y_pred = model.predict(x)
	# print(y_pred[0])

	# Save
	with open('model{}.json'.format(tag), 'w') as outfile:
		outfile.write(model.to_json()) 
	if hist:
		with open('hist{}.csv'.format(tag), 'w') as outfile:
			outfile.write(','.join(hist.history.keys()) + '\n')
			for i in range(len(hist.history.values()[0])):
				outfile.write(','.join([str(x[i]) for x in hist.history.values()]) + '\n')
	else:
		print('No history to save')

