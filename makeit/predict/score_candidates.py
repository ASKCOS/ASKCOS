# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import numpy as np
import os
import sys
import argparse
from keras import backend as K 
from keras.models import Sequential, model_from_json
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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt    # for visualization


def build(N_e = 2048, N_h1 = 100, N_h2 = 50, N_h3 = 0, N_c = 500, inner_act = 'tanh',
		dr = 0.5, l2v = 0.01, lr = 0.0003):
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
		Dense(N_h1, activation = inner_act, W_regularizer = l2(l2v)), 
		input_shape = (N_c, N_e)
	))
	model.add(Dropout(dr))

	# Time distributed dense
	if N_h2 > 0:
		model.add(TimeDistributed(
			Dense(N_h2, activation = inner_act, W_regularizer = l2(l2v))
		))
		model.add(Dropout(dr))

	# Time distributed dense
	if N_h3 > 0:
		model.add(TimeDistributed(
			Dense(N_h3, activation = inner_act, W_regularizer = l2(l2v))
		))
		model.add(Dropout(dr))

	# Time distributed dense
	model.add(TimeDistributed(
		Dense(1, activation = 'linear', W_regularizer = l2(l2v))
	))

	# Flatten to 1 x N_c
	model.add(Flatten())

	# Softmax
	model.add(Activation('softmax'))

	# Now compile
	sgd = SGD(lr = lr, decay = 1e-4, momentum = 0.9)
	# model.compile(loss = 'binary_crossentropy', optimizer = sgd, 
	# 	metrics = ['binary_accuracy'])
	model.compile(loss = 'categorical_crossentropy', optimizer = sgd, 
		metrics = ['accuracy'])

	return model

def train(model, x_files, y_files, tag = '', split_ratio = 0.8):

	hist_fid = open(os.path.join(FROOT, 'hist{}.csv'.format(tag)), 'a')
	hist_fid.write('epoch,filenum,loss,val_loss,acc,val_acc\n')
	try:
		for epoch in range(nb_epoch):
			print('>>> Epoch {}/{} <<<'.format(epoch + 1, nb_epoch))
			# Get each set of data
			for fnum in range(len(x_files)):
				print('Running file set {} for training/validation'.format(fnum))
				x = get_x_data(x_files[fnum])
				y = get_y_data(y_files[fnum])

				hist = model.fit(x, y, 
					nb_epoch = 1, 
					batch_size = batch_size, 
					validation_split = (1 - split_ratio),
					verbose = False)

				hist_fid.write('{},{},{},{},{},{}\n'.format(
					epoch + 1, fnum, hist.history['loss'][0], hist.history['val_loss'][0],
					hist.history['acc'][0], hist.history['val_acc'][0]
				))
				print('loss: {}, val_loss: {}, acc: {}, val_acc: {}'.format(
					hist.history['loss'][0], hist.history['val_loss'][0],
					hist.history['acc'][0], hist.history['val_acc'][0]
				))
	except KeyboardInterrupt:
		print('Stopped training early!')
		hist_fid.close()
		return None

	hist_fid.close()

	return True

def get_x_data(fpath):
	with open(fpath, 'rb') as infile:
		x = pickle.load(infile)
		x = [np.array(coo_matrix((a[0], (a[1], a[2])), shape = a[3]).todense()) for a in x]
		x = np.transpose(np.array(x), (0, 2, 1))
		if USE_EXPLICIT_DIFF:
			h = x.shape[2] / 2 # halfway in concatenated FP
			return np.concatenate((x, x[:,:,h:] - x[:,:,:h]), axis = 2)
		return x

def get_y_data(fpath):
	with open(fpath, 'rb') as infile:
		y = pickle.load(infile)
		y = np.array(coo_matrix((y[0], (y[1], y[2])), shape = y[3]).todense())
		return y

def get_z_data(fpath):
	with open(fpath, 'rb') as infile:
		z = pickle.load(infile)
		return z

def pred_histogram(model, x_files, y_files, z_files, tag = '', split_ratio = 0.8):
	'''
	Given a trained model and a list of samples, this function creates a 
	histogram of the pseudo-probabilities assigned to the true reaction 
	examples (i.e., the correct product)
	'''

	fid = open(os.path.join(FROOT, 'probs{}.dat'.format(tag)), 'w')

	# Calculate probabilities
	
	train_preds = []; val_preds = []
	train_corr = 0; val_corr = 0;
	dataset = ''; trueprob = 0; trueprobs = []

	for fnum in range(len(x_files)):
		print('Testing file number {}'.format(fnum))
		# Data must be pre-padded
		x = get_x_data(x_files[fnum])
		y = get_y_data(y_files[fnum])
		z = get_z_data(z_files[fnum])

		preds = model.predict(x, batch_size = batch_size)
		trueprobs = []

		for i in range(preds.shape[0]): 
			pred = preds[i, :] # Iterate through each sample
			trueprob = pred[y[i,:]] # prob assigned to true outcome
			trueprobs.append(trueprob)
			if i < int(split_ratio * preds.shape[0]):
				dataset = 'train'
				train_preds.append(trueprob)
				if np.argmax(pred) == np.argmax(y[i,:]):
					train_corr += 1
			else:
				dataset = 'val'
				val_preds.append(trueprob)
				if np.argmax(pred) == np.argmax(y[i,:]):
					val_corr += 1
			fid.write('{}\t{}\t{}\n'.format(z[i], dataset, trueprobs[i][0]))
	fid.close()
	
	train_acc = train_corr / float(len(train_preds))
	val_acc = val_corr / float(len(val_preds))

	train_preds = np.array(train_preds)
	val_preds = np.array(val_preds)


	def histogram(array, title, label, acc):
		acc = int(acc * 1000)/1000. # round3
		try:
			# Visualize in histogram
			weights = np.ones_like(array) / len(array)
			plt.clf()
			n, bins, patches = plt.hist(array, np.arange(0, 1, 0.02), facecolor = 'blue', alpha = 0.5, weights = weights)
			plt.xlabel('Assigned probability to true product')
			plt.ylabel('Normalized frequency')
			plt.title('Histogram of pseudo-probabilities - {} (N={},acc={})'.format(title, len(array), acc))
			plt.axis([0, 1, 0, 1])
			plt.grid(True)
			plt.savefig(os.path.join(FROOT, label), bbox_inches = 'tight')
		except:
			pass

	histogram(train_preds, 'TRAIN', 'training{}.png'.format(tag), train_acc)
	histogram(val_preds, 'VAL', 'validation{}.png'.format(tag), val_acc)

def visualize_weights(model, tag):
	'''Given a trained model, this function visualizes the weights in each Dense layer'''
	
	def vec_to_png(vec, label, title):
		try:
			print(vec.shape)
			fig = plt.figure(figsize=(4,4))
			plt.pcolor(vec)#, vmin = 0, vmax = 1, cmap = plt.get_cmap('Greens'))
			plt.title('{} {}'.format(title, vec.shape))
			cbar = plt.colorbar()
			#plt.gca().yaxis.set_visible(False)
			#plt.gca().xaxis.set_visible(False)
			plt.xlim([0, vec.shape[1]])
			plt.ylim([0, vec.shape[0]])
			#plt.subplots_adjust(left = 0, right = 1, top = 0.4, bottom = 0)
			plt.savefig(os.path.join(FROOT, label) + '.png', bbox_inches = 'tight')
			plt.close(fig)
			plt.clf()
		except:
			pass

	for i, layer in enumerate(model.layers):
		try:
			weights = layer.layer.W.eval()
			vec_to_png(weights, 'weights{}(layer{})'.format(tag, i), 'Layer {}'.format(i))
		except Exception as e:
			print('Skipping layer {}'.format(i))
			print(e)


if __name__ == '__main__':

	FROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'output')
	
	parser = argparse.ArgumentParser()
	parser.add_argument('--nb_epoch', type = int, default = 100,
						help = 'Number of epochs to train for')
	parser.add_argument('--batch_size', type = int, default = 40,
						help = 'Batch size')
	parser.add_argument('--Nh1', type = int, default = 40,
						help = 'Number of hidden nodes in first layer')
	parser.add_argument('--Nh2', type = int, default = 0,
						help = 'Number of hidden nodes in second layer')
	parser.add_argument('--Nh3', type = int, default = 0,
						help = 'Number of hidden nodes in third layer')
	parser.add_argument('--tag', type = str, default = int(time.time()),
						help = 'Tag for this model')
	parser.add_argument('--retrain', type = bool, default = False,
		                help = 'Retrain with loaded weights')
	parser.add_argument('--test', type = bool, default = False,
						help = 'Test model only')
	parser.add_argument('--l2', type = float, default = 0.01,
						help = 'l2 regularization parameter for each Dense layer, default 0.01')
	parser.add_argument('data', type = str,
		                help = 'Data folder with x*, y*, and z* data files')
	parser.add_argument('--lr', type = float, default = 0.0003, 
						help = 'Learning rate, default 0.0003')
	parser.add_argument('--dr', type = float, default = 0.5,
						help = 'Dropout rate, default 0.5')
	parser.add_argument('--visualize', type = bool, default = False,
		                help = 'Whether or not to visualize weights ONLY, default False')
	parser.add_argument('--explicitdiff', type = bool, default = False,
						help = 'Include explicit fingerprint difference in X, default False')
	args = parser.parse_args()

	x_files = sorted([os.path.join(args.data, dfile) \
					for dfile in os.listdir(args.data) \
					if dfile[0] == 'x'])
	y_files = sorted([os.path.join(args.data, dfile) \
					for dfile in os.listdir(args.data) \
					if dfile[0] == 'y'])
	z_files = sorted([os.path.join(args.data, dfile) \
					for dfile in os.listdir(args.data) \
					if dfile[0] == 'z'])
	print(x_files)
	print(y_files)
	print(z_files)

	nb_epoch = int(args.nb_epoch)
	batch_size = int(args.batch_size)
	N_h1 = int(args.Nh1)
	N_h2 = int(args.Nh2)
	N_h3 = int(args.Nh3)
	l2v = float(args.l2)
	lr = float(args.lr)
	dr = float(args.dr)
	N_c = 500
	N_e = 2048

	tag = args.tag

	USE_EXPLICIT_DIFF = False
	if bool(args.explicitdiff):
		N_e = 3 * N_e / 2
		USE_EXPLICIT_DIFF = True

	if bool(args.retrain):
		model = model_from_json(open(os.path.join(FROOT, 'model{}.json'.format(tag))).read())
		model.compile(loss = 'categorical_crossentropy', 
			optimizer = SGD(lr = lr, decay = 1e-4, momentum = 0.9),
			metrics = ['accuracy']
		)
		model.load_weights(os.path.join(FROOT, 'weights{}.h5'.format(tag)))
	else:
		model = build(N_e = N_e, N_c = N_c, N_h1 = N_h1, N_h2 = N_h2, N_h3 = N_h3, l2v = l2v, lr = lr, dr = dr)
	
	if bool(args.test):
		pred_histogram(model, x_files, y_files, z_files, tag = tag, split_ratio = 0.8)
		quit(1)
	elif bool(args.visualize):
		visualize_weights(model, tag)
		quit(1)

	hist = train(model, x_files, y_files, tag = tag, split_ratio = 0.8)
	model.save_weights(os.path.join(FROOT, 'weights{}.h5'.format(tag)))

	# Save
	with open(os.path.join(FROOT, 'model{}.json'.format(tag)), 'w') as outfile:
		outfile.write(model.to_json()) 

	pred_histogram(model, x_files, y_files, z_files, tag = tag, split_ratio = 0.8)
