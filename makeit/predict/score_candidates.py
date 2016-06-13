# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import numpy as np
import os
from keras import backend as K 
from keras.models import Sequential
from keras.layers import Dense, Activation
from keras.layers.core import Flatten, Permute, Reshape, Dropout
from keras.optimizers import *
from keras.layers.convolutional import Convolution1D, Convolution2D
from makeit.predict.preprocess_candidates import *
import theano.tensor as T
from scipy.sparse import coo_matrix
import cPickle as pickle

def build(N_e = 2048, N_h1 = 50, N_h2 = 20, N_c = 10000, inner_act = 'relu',
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

	# Change shape to 1 x 2048 x 10000
	model.add(Reshape((1, N_e, N_c), input_shape = (N_e, N_c)))

	# 2048-long convolution (N_h1 filters)
	model.add(Convolution2D(N_h1, N_e, 1, input_shape = (1, N_e, N_c) ))

	# Dropout
	model.add(Dropout(dr))

	# Activation
	model.add(Activation(inner_act))

	# Reshape again to 1 x 100 x 10000
	model.add(Permute((2, 1, 3), input_shape = (N_h1, 1, N_c)))

	# Another convolution (N_h2 filters)
	model.add(Convolution2D(N_h2, N_h1, 1, input_shape = (1, N_h1, N_c) ))

	# Dropout
	model.add(Dropout(dr))

	# Activation
	model.add(Activation(inner_act))

	# Reshape again to 1 x 50 x 10000
	model.add(Permute((2, 1, 3), input_shape = (N_h2, 1, N_c)))

	# Another convolution (1 filter)
	model.add(Convolution2D(1, N_h2, 1, input_shape = (1, N_h2, N_c) ))

	# Dropout
	model.add(Dropout(dr))

	# Flatten to 1 x N_c
	model.add(Flatten(input_shape = (1, 1, N_c)))

	# Add on a softmax to to "select" the most likely candidate
	model.add(Activation('sigmoid'))

	# Now compile
	sgd = SGD(lr = 0.09, decay = 1e-6, momentum = 0.9)
	model.compile(loss = 'binary_crossentropy', optimizer = sgd, 
		metrics = ['binary_accuracy'])

	return model

def train(model, x, y):

	nb_epoch = 1000
	batch_size = 15

	hist = model.fit(x, y, 
		nb_epoch = nb_epoch, 
		batch_size = batch_size, 
		validation_split = 0.2,
		verbose = 1)

	return hist


if __name__ == '__main__':

	# Data must be pre-padded
	with open('x_coo.dat', 'rb') as infile:
		x = pickle.load(infile)
		x = [np.array(coo_matrix((z[0], (z[1], z[2])), shape = z[3]).todense()) for z in x]
		x = np.array(x)
	with open('y_coo.dat', 'rb') as infile:
		y = pickle.load(infile)
		y = np.array(coo_matrix((y[0], (y[1], y[2])), shape = y[3]).todense())

	# Build and train
	N_c = len(y[0,:])
	print('N_c: {}'.format(N_c))
	print('Mean number of true predictions: {}'.format(np.mean(np.sum(y, axis = 1))))

	model = build(N_c = N_c)
	hist = train(model, x, y)

	y_pred = model.predict(x)
	print(y_pred[0])

	with open('hist.csv', 'w') as outfile:
		outfile.write(','.join(hist.history.keys()) + '\n')
		for i in range(len(hist.history.values()[0])):
			outfile.write(','.join([str(x[i]) for x in hist.history.values()]) + '\n')
