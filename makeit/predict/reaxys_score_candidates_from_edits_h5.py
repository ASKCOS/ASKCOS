# Import relevant packages
from __future__ import print_function
import time
from global_config import USE_STEREOCHEMISTRY
import numpy as np
import os
import sys
import argparse
import h5py # needed for save_weights, fails otherwise
from keras import backend as K 
import theano
from keras.models import Sequential, Model, model_from_json
from keras.layers import Dense, Activation, Input
from keras.layers.core import Flatten, Permute, Reshape, Dropout, Lambda, RepeatVector
from keras.layers.wrappers import TimeDistributed
from keras.engine.topology import Merge, merge
from keras.optimizers import *
from keras.layers.convolutional import Convolution1D, Convolution2D
from keras.regularizers import l2
from keras.utils.np_utils import to_categorical
no_printing = False
try:
	from keras.utils.visualize_util import plot
except:
	no_printing = True
from makeit.embedding.descriptors import edits_to_vectors, oneHotVector # for testing
import rdkit.Chem as Chem
import theano.tensor as T
from scipy.sparse import coo_matrix
import cPickle as pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt    # for visualization
import scipy.stats as ss
import itertools


def build(F_atom = 1, F_bond = 1, N_e = 5, N_h1 = 100, N_h2 = 50, N_h3 = 0, N_c = 100, inner_act = 'tanh',
		l2v = 0.01, lr = 0.0003, N_hf = 20, context_weight = 150.0, enhancement_weight = 0.1):
	'''
	Builds the feed forward model.

	N_e:  maximum number of edits of each type
	N_h1: number of hidden nodes in first layer
	N_h2: number of hidden nodes in second layer
	N_c:  maximum number of candidates for each sample (padded length)
	inner_act: activation function 
	'''

	h_lost = Input(shape = (N_c, N_e, F_atom), name = "H_lost")
	h_gain = Input(shape = (N_c, N_e, F_atom), name = "H_gain")
	bond_lost = Input(shape = (N_c, N_e, F_bond), name = "bond_lost")
	bond_gain = Input(shape = (N_c, N_e, F_bond), name = "bond_gain")
	reagents = Input(shape = (256,), name = "reagent FP") # TODO: remove hard-coded length
	solvent = Input(shape = (6,), name = "solvent descriptors c,e,s,a,b,v")
	temp = Input(shape = (1,), name = "temperature [C]")

	h_lost_r = Reshape((N_c*N_e, F_atom), name = "flatten H_lost")(h_lost)
	h_gain_r = Reshape((N_c*N_e, F_atom), name = "flatten H_gain")(h_gain)
	bond_lost_r = Reshape((N_c*N_e, F_bond), name = "flatten bond_lost")(bond_lost)
	bond_gain_r = Reshape((N_c*N_e, F_bond), name = "flatten bond_gain")(bond_gain)

	h_lost_h1 = TimeDistributed(Dense(N_h1, activation = inner_act, W_regularizer = l2(l2v)), name = "embed H_lost 1")(h_lost_r)
	h_gain_h1 = TimeDistributed(Dense(N_h1, activation = inner_act, W_regularizer = l2(l2v)), name = "embed H_gain 1")(h_gain_r)
	bond_lost_h1 = TimeDistributed(Dense(N_h1, activation = inner_act, W_regularizer = l2(l2v)), name = "embed bond_lost 1")(bond_lost_r)
	bond_gain_h1 = TimeDistributed(Dense(N_h1, activation = inner_act, W_regularizer = l2(l2v)), name = "embed bond_gain 1")(bond_gain_r)
	N_h = N_h1

	if N_h2 > 0:
		h_lost_h2 = TimeDistributed(Dense(N_h2, activation = inner_act, W_regularizer = l2(l2v)), name = "embed H_lost 2")(h_lost_h1)
		h_gain_h2 = TimeDistributed(Dense(N_h2, activation = inner_act, W_regularizer = l2(l2v)), name = "embed H_gain 2")(h_gain_h1)
		bond_lost_h2 = TimeDistributed(Dense(N_h2, activation = inner_act, W_regularizer = l2(l2v)), name = "embed bond_lost 2")(bond_lost_h1)
		bond_gain_h2 = TimeDistributed(Dense(N_h2, activation = inner_act, W_regularizer = l2(l2v)), name = "embed bond_gain 2")(bond_gain_h1)
		N_h = N_h2

		if N_h3 > 0:
			h_lost_h = TimeDistributed(Dense(N_h3, activation = inner_act, W_regularizer = l2(l2v)), name = "embed H_lost 3")(h_lost_h2)
			h_gain_h = TimeDistributed(Dense(N_h3, activation = inner_act, W_regularizer = l2(l2v)), name = "embed H_gain 3")(h_gain_h2)
			bond_lost_h = TimeDistributed(Dense(N_h3, activation = inner_act, W_regularizer = l2(l2v)), name = "embed bond_lost 3")(bond_lost_h2)
			bond_gain_h = TimeDistributed(Dense(N_h3, activation = inner_act, W_regularizer = l2(l2v)), name = "embed bond_gain 3")(bond_gain_h2)
			N_h = N_h3

		else:
			h_lost_h = h_lost_h2
			h_gain_h = h_gain_h2
			bond_lost_h = bond_lost_h2
			bond_gain_h = bond_gain_h2

	else:
		h_lost_h = h_lost_h1
		h_gain_h = h_gain_h1
		bond_lost_h = bond_lost_h1
		bond_gain_h = bond_gain_h1

	h_lost_r2 = Reshape((N_c, N_e, N_h), name = "expand H_lost edits")(h_lost_h)
	h_gain_r2 = Reshape((N_c, N_e, N_h), name = "expand H_gain edits")(h_gain_h)
	bond_lost_r2 = Reshape((N_c, N_e, N_h), name = "expand bond_lost edits")(bond_lost_h)
	bond_gain_r2 = Reshape((N_c, N_e, N_h), name = "expand bond_gain edits")(bond_gain_h)

	h_lost_sum = Lambda(lambda x: K.sum(x, axis = 2), output_shape = (N_c, N_h), name = "sum H_lost")(h_lost_r2)
	h_gain_sum = Lambda(lambda x: K.sum(x, axis = 2), output_shape = (N_c, N_h), name = "sum H_gain")(h_gain_r2)
	bond_lost_sum = Lambda(lambda x: K.sum(x, axis = 2), output_shape = (N_c, N_h), name = "sum bond_lost")(bond_lost_r2)
	bond_gain_sum = Lambda(lambda x: K.sum(x, axis = 2), output_shape = (N_c, N_h), name = "sum bond_gain")(bond_gain_r2)

	net_sum = merge([h_lost_sum, h_gain_sum, bond_lost_sum, bond_gain_sum], mode = 'sum', name = "sum across edits")

	feature_to_feature = Dense(N_hf, activation = inner_act, W_regularizer = l2(l2v))
	net_sum_h = TimeDistributed(feature_to_feature, name = "reaction embedding post-sum")(net_sum)

	# Take reagents -> intermediate representation -> cosine similarity to enhance reaction
	reagents_h = Dense(N_hf, activation = 'tanh', W_regularizer = l2(l2v), name = "reagent fingerprint to features")(reagents)
	reagents_h_rpt = RepeatVector(N_c, name = "broadcast reagent vector")(reagents_h)

	# Dot product between reagents and net_sum_h gives enhancement factor
	enhancement_mul = merge([net_sum_h, reagents_h_rpt], mode = 'mul', name = "multiply reaction with reagents [dot 1/2]")
	enhancement = Lambda(lambda x: K.sum(x, axis = -1), output_shape = (N_c, 1), name = "sum reaction with reagents [dot 2/2]")(enhancement_mul)
	enhancement_r = Reshape((N_c, 1), name = "shape check 1")(enhancement) # needs to be explicit for some reason

	# Converge to G0, C[not real], E, S, A, B, V, and K
	feature_to_params = Dense(8, activation = 'linear', W_regularizer = l2(l2v))
	params = TimeDistributed(feature_to_params, name = "features to K,G0,C,E,S,A,B,V")(net_sum_h)

	# Concatenate enhancement and solvents
	solvent_rpt = RepeatVector(N_c, name = "broadcast solvent vector")(solvent)
	temp_rpt = RepeatVector(N_c, name = "broadcast temperature")(temp)
	params_enhancement = merge([params, enhancement_r, solvent_rpt, temp_rpt], mode = 'concat', name = "concatenate context")

	# # Calculate using thermo-ish
	# # K * exp(- (G0 + delG_solv) / T + enhancement)
	# unscaled_score = Lambda(
	# 	lambda x: x[:, :, 0] * K.exp(- (x[:, :, 1] + K.sum(x[:, :, 2:8] * x[:, :, 8:14], axis = -1)) / (x[:, :, 15] + 273.15) + x[:, :, 8]),
	# 	output_shape = lambda x: (None, N_c,),
	# 	name = "propensity = K * exp(- (G0 + cC + eE + ... + vV) / T + enh.)"
	# )(params_enhancement)

	#### NON-EXPONENTIAL VERSION
	unscaled_score = Lambda(
		lambda x: x[:, :, 0] - context_weight * (x[:, :, 1] + K.sum(x[:, :, 2:8] * x[:, :, 9:15], axis = -1)) / (x[:, :, 15] + 273.15) + enhancement_weight * x[:, :, 8],
		output_shape = lambda x: (None, N_c,),
		name = "propensity = logK - (G0 + cC + eE + ... + vV) / T + enh."
	)(params_enhancement)

	unscaled_score_r = Reshape((N_c,), name = "shape check 2")(unscaled_score)

	score = Activation('softmax', name = "scores to probs")(unscaled_score_r)
	#score = unscaled_score_r

	model = Model(input = [h_lost, h_gain, bond_lost, bond_gain, reagents, solvent, temp], 
		output = [score])

	model.summary()

	# Now compile
	sgd = SGD(lr = lr, decay = 1e-4, momentum = 0.9)
	adam = Adam(lr = lr)

	model.compile(loss = 'categorical_crossentropy', optimizer = adam, 
		metrics = ['accuracy'])

	return model

def data_generator(h5f, start_at = 0, end_at = None, batch_size = 20):
	'''This function generates batches of data from the
	h5f file since all the data can't fit in memory.

	The starting and ending indices are specified explicitly so the
	same function can be used for validation data as well'''

	inputs = [
		'x_h_lost', 
		'x_h_gain',
		'x_bond_lost', 
		'x_bond_gain', 
		'reagent', 
		'solvent', 
		'T'
	]
	outputs = [
		'reaction_true_onehot'
	]

	if type(end_at) == type(None):
		end_at = len(h5f['x_h_lost'])

	# Keep returning forever and ever
	while True:
		for startIndex in range(start_at, end_at, batch_size):
			endIndex = min(startIndex + batch_size, end_at)
			# print('Batch {} to {}'.format(startIndex, endIndex))
			# yield (x, y) as tuple, but each one is a list
			yield (
				[h5f[dset][startIndex:endIndex] for dset in inputs],
				[h5f[dset][startIndex:endIndex] for dset in outputs],
			)


def train(model, h5f, tag = '', split_ratio = 0.8):

	# hist_fid = open(os.path.join(FROOT, 'hist{}.csv'.format(tag)), 'a')
	# hist_fid.write('epoch,filenum,loss,val_loss,acc,val_acc\n')
	
	N_samples =  len(h5f['x_h_lost'])
	N_train = int(N_samples * split_ratio)
	print('Total number of samples: {}'.format(N_samples))
	print('Training on {}% - {}'.format(split_ratio*100, N_train))

	try:

		hist = model.fit_generator(
			data_generator(h5f, 0, N_train, batch_size), 
			samples_per_epoch = N_train,
			nb_epoch = nb_epoch, 
			validation_data = data_generator(h5f, N_train, N_samples, batch_size),
			nb_val_samples = N_samples - N_train,
			#pickle_safe = True
			verbose = 1,
		)

		# hist_fid.write('{},{},{},{},{},{}\n'.format(
		# 	epoch + 1, fnum, hist.history['loss'][0], hist.history['val_loss'][0],
		# 	hist.history['acc'][0], hist.history['val_acc'][0]
		# ))
		# print('loss: {}, val_loss: {}, acc: {}, val_acc: {}'.format(
		#	hist.history['loss'][0], hist.history['val_loss'][0],
		# 	hist.history['acc'][0], hist.history['val_acc'][0]
		#))
		# 	average_acc += hist.history['acc'][0]
		# 	average_val_loss += hist.history['val_loss'][0]
		# 	average_val_acc += hist.history['val_acc'][0]
		# print('    loss:     {:8.4f}, acc:     {:5.4f}'.format(average_loss/len(x_files), average_acc/len(x_files)))
		# print('    val_loss: {:8.4f}, val_acc: {:5.4f}'.format(average_val_loss/len(x_files), average_val_acc/len(x_files)))
		model.save_weights(os.path.join(FROOT, 'weights{}.h5'.format(tag)), overwrite = True)
	except KeyboardInterrupt:
		print('Stopped training early!')
		# hist_fid.close()
		return None

	# hist_fid.close()

	return True


def pred_histogram(model, x_files, xc_files, y_files, z_files, tag = '', split_ratio = 0.8):
	'''
	Given a trained model and a list of samples, this function creates a 
	histogram of the pseudo-probabilities assigned to the true reaction 
	examples (i.e., the correct product)
	'''

	fid = open(os.path.join(FROOT, 'probs{}.dat'.format(tag)), 'w')
	fid.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
		'reaction_smiles', 'train/val', 
		'true_edit', 'prob_true_edit', 
		'predicted_edit(or no. 2)', 'prob_predicted_edit(or no. 2)',
		'rank_true_edit'
	))

	# Calculate probabilities
	
	train_preds = []; val_preds = []
	train_corr = 0; val_corr = 0;
	dataset = ''; trueprob = 0; trueprobs = []

	for fnum in range(len(x_files)):
		print('Testing file number {}'.format(fnum))
		# Data must be pre-padded
		z = get_z_data(z_files[fnum])
		x = list(get_x_data(x_files[fnum], z))
		xc = list(get_xc_data(xc_files[fnum]))
		y = get_y_data(y_files[fnum])

		# Getting unprocessed x data
		with open(x_files[fnum], 'rb') as infile:
			all_edits = pickle.load(infile)

		preds = model.predict(x + xc, batch_size = 20)
		trueprobs = []

		for i in range(preds.shape[0]): 
			rank_true_edit = 0
			edits = all_edits[i]
			pred = preds[i, :] # Iterate through each sample
			trueprob = pred[y[i,:] != 0][0] # prob assigned to true outcome
			trueprobs.append(trueprob)
			rank_true_edit = 1 + len(pred) - (ss.rankdata(pred))[np.argmax(y[i,:])]
			
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

			if rank_true_edit != 1:
				# record highest probability
				most_likely_edit_i = np.argmax(pred)
				most_likely_prob = np.max(pred)
			else:
				# record number two prediction
				most_likely_edit_i = np.argmax(pred[pred != np.max(pred)])
				most_likely_prob = np.max(pred[pred != np.max(pred)])

			if most_likely_edit_i >= len(edits): # no reaction?
				most_likely_edit = 'no_reaction'
			else:
				most_likely_edit = edits[most_likely_edit_i]
			fid.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
				z[i], dataset, 
				edits[np.argmax(y[i,:])], trueprob, 
				most_likely_edit, most_likely_prob,
				rank_true_edit
			))
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
			n, bins, patches = plt.hist(array, np.arange(0, 1.02, 0.02), facecolor = 'blue', alpha = 0.5, weights = weights)
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


if __name__ == '__main__':

	FROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'output')
	
	parser = argparse.ArgumentParser()
	parser.add_argument('--nb_epoch', type = int, default = 100,
						help = 'Number of epochs to train for, default 100')
	parser.add_argument('--batch_size', type = int, default = 20,
						help = 'Batch size, default 20')
	parser.add_argument('--Nh1', type = int, default = 40,
						help = 'Number of hidden nodes in first layer, default 40')
	parser.add_argument('--Nh2', type = int, default = 0,
						help = 'Number of hidden nodes in second layer, default 0')
	parser.add_argument('--Nh3', type = int, default = 0,
						help = 'Number of hidden nodes in third layer, ' + 
								'immediately before summing, default 0')
	parser.add_argument('--Nhf', type = int, default = 20,
						help = 'Number of hidden nodes in layer between summing ' +
								'and final score, default 20')
	parser.add_argument('--tag', type = str, default = str(time.time()),
						help = 'Tag for this model')
	parser.add_argument('--retrain', type = bool, default = False,
		                help = 'Retrain with loaded weights, default False')
	parser.add_argument('--test', type = bool, default = False,
						help = 'Test model only, default False')
	parser.add_argument('--l2', type = float, default = 0.01,
						help = 'l2 regularization parameter for each Dense layer, default 0.01')
	parser.add_argument('data', type = str,
		                help = 'Data file')
	parser.add_argument('--lr', type = float, default = 0.01, 
						help = 'Learning rate, default 0.01')
	# parser.add_argument('--dr', type = float, default = 0.5,
	# 					help = 'Dropout rate, default 0.5')
	parser.add_argument('--fold', type = int, default = 5, 
						help = 'Which fold of the 5-fold CV is this? Defaults 5')
	parser.add_argument('--visualize', type = bool, default = False,
				help = 'Whether or not to visualize weights ONLY, default False')
	parser.add_argument('--context_weight', type = float, default = 100.0,
					help = 'Weight assigned to contextual effects, default 100.0')
	parser.add_argument('--enhancement_weight', type = float, default = 0.1,
			help = 'Weight assigned to enhancement factor, default 0.1')
	parser.add_argument('--Nc', type = int, default = 500,
			help = 'Number of candidates per example, default 500')

	args = parser.parse_args()

	mol = Chem.MolFromSmiles('[C:1][C:2]')
	(a, _, b, _) = edits_to_vectors((['1'],[],[('1','2',1.0)],[]), mol)

	F_atom = len(a[0])
	F_bond = len(b[0])

	nb_epoch = int(args.nb_epoch)
	batch_size = int(args.batch_size)
	N_h1 = int(args.Nh1)
	N_h2 = int(args.Nh2)
	N_h3 = int(args.Nh3)
	N_hf = int(args.Nhf)
	l2v = float(args.l2)
	lr = float(args.lr)
	N_c = int(args.Nc) # number of candidate edit sets
	N_e = 5 # maximum number of edits per class
	context_weight = float(args.context_weight)
	enhancement_weight = float(args.enhancement_weight)
	data_file = args.data

	THIS_FOLD_OUT_OF_FIVE = int(args.fold)
	tag = args.tag + ' fold{}'.format(args.fold)

	if bool(args.retrain):
		print('Reloading from file')
		rebuild = raw_input('Do you want to rebuild from scratch instead of loading from file? [n/y] ')
		if rebuild == 'y':
			model = build(F_atom = F_atom, F_bond = F_bond, N_e = N_e, N_c = N_c, N_h1 = N_h1, 
				N_h2 = N_h2, N_h3 = N_h3, N_hf = N_hf, l2v = l2v, lr = lr, 
				context_weight = context_weight, enhancement_weight = enhancement_weight)
		else:
			model = model_from_json(open(os.path.join(FROOT, 'model{}.json'.format(tag))).read())
			model.compile(loss = 'categorical_crossentropy', 
				optimizer = Adam(lr = lr),
				metrics = ['accuracy'])
			
		model.load_weights(os.path.join(FROOT, 'weights{}.h5'.format(tag)))
	else:
		model = build(F_atom = F_atom, F_bond = F_bond, N_e = N_e, N_c = N_c, N_h1 = N_h1, N_h2 = N_h2, N_h3 = N_h3, N_hf = N_hf, l2v = l2v, lr = lr, context_weight = context_weight)
		try:
			with open(os.path.join(FROOT, 'model{}.json'.format(tag)), 'w') as outfile:
				outfile.write(model.to_json())
		except:
			print('could not write model to json')

	if bool(args.test):
		pred_histogram(model, x_files, xc_files, y_files, z_files, tag = tag, split_ratio = 0.8)
		quit(1)
	elif bool(args.visualize):
		visualize_weights(model, tag)
		quit(1)

	h5f = h5py.File(data_file, 'r')
	hist = train(model, h5f, tag = tag, split_ratio = 0.8)
	model.save_weights(os.path.join(FROOT, 'weights{}.h5'.format(tag)), overwrite = True) 

	#pred_histogram(model, x_files, xc_files, y_files, z_files, tag = tag, split_ratio = 0.8)
	#visualize_weights(model, tag)
	h5f.close()