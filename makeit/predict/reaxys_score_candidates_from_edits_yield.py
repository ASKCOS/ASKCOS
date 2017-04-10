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

def msle_of_true(y_true, y_pred):
	'''Custom loss function that uses the mean squared log error in predicted
	yield, assuming that y_pred are unscaled predictions with the true 
	outcome in the first index.'''
	return K.square(K.log(K.clip(y_pred[:, 0:1], K.epsilon(), 1.0)) - K.log(K.clip(y_true, K.epsilon(), 1.0)))

def mse_of_true(y_true, y_pred):
	'''Custom loss function that uses the mean squared error in predicted
	yield, assuming that y_pred are unscaled predictions with the true 
	outcome in the first index.'''
	return K.square(y_pred[:, 0:1] - y_true)


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

	scaled_score = Activation(lambda x: K.exp(x / 5.0 - 1.0), name = 'exponential activation')(unscaled_score_r)

	# Do not scale score with softmax (which would force 100% conversion)
	# Scale linearly

	yields = Lambda(
		lambda x: x / K.tile(K.maximum(1.0, K.sum(x, axis = -1, keepdims = True)), (1, N_c)),
	name = "scale if sum(yields)>1")(scaled_score)

	model = Model(input = [h_lost, h_gain, bond_lost, bond_gain, reagents, solvent, temp], 
		output = [yields])

	model.summary()
	# if no_printing:
	# 	print('Could not print')
	# else:
	# 	plot(model, to_file = 'model.png', show_shapes = True)


	# Now compile
	adam = Adam(lr = lr)

	model.compile(loss = mse_of_true, optimizer = adam)

	return model

def train(model, x_files, xc_files, y_files, z_files, tag = '', split_ratio = 0.8):

	hist_fid = open(os.path.join(FROOT, 'hist{}.csv'.format(tag)), 'a')
	hist_fid.write('epoch,filenum,loss,val_loss\n')
	try:
		for epoch in range(nb_epoch):
			print('>>> EPOCH {}/{} <<<'.format(epoch + 1, nb_epoch))
			average_loss = 0;
			average_val_loss = 0;
			# Shuffle order that files are read
			fnums = range(len(x_files))
			np.random.shuffle(fnums)
			for zz, fnum in enumerate(fnums): # get each set of data
				print('    running file {}/{} (no. {})'.format(fnum+1, len(x_files), zz+1))
				
				if len(x_files) > 1 or epoch == 0: # get new data
					z = get_z_data(z_files[fnum])
					x = list(get_x_data(x_files[fnum], z))
					xc = list(get_xc_data(xc_files[fnum]))
					y = get_y_data(y_files[fnum])
					z = rearrange_for_5fold_cv(z)
					
				# print(x[0].shape)
				# print(x[1].shape)
				# print(x[2].shape)
				# print(x[3].shape)
				# print(xc[0].shape)
				# print(xc[1].shape)
				# print(xc[2].shape)
				# print('y: {}'.format(y.shape))

				# ys = model.predict(x + xc)
				# print(ys)
				# print(ys.shape)
				# mse = msle_of_true(y, ys)
				# print('MSE: {}'.format(mse.eval()))
				# raw_input('Predicted, pause...')

				hist = model.fit(x + xc, y, 
					nb_epoch = 1, 
					batch_size = batch_size, 
					validation_split = (1 - split_ratio),
					verbose = 0)

				hist_fid.write('{},{},{},{}\n'.format(
					epoch + 1, fnum, hist.history['loss'][0], hist.history['val_loss'][0],
				))
				# print('loss: {}, val_loss: {}, acc: {}, val_acc: {}'.format(
				#	hist.history['loss'][0], hist.history['val_loss'][0],
				# 	hist.history['acc'][0], hist.history['val_acc'][0]
				#))
	                        print('This train loss: {}'.format(hist.history['loss'][0]))
                                print('This val loss:   {}'.format(hist.history['val_loss'][0]))
                                average_loss += hist.history['loss'][0]
				average_val_loss += hist.history['val_loss'][0]
			print('    loss:     {:8.4f}'.format(average_loss/len(x_files)))
			print('    val_loss: {:8.4f}'.format(average_val_loss/len(x_files)))
			model.save_weights(os.path.join(FROOT, 'weights{}.h5'.format(tag)), overwrite = True)
	except KeyboardInterrupt:
		print('Stopped training early!')
		hist_fid.close()
		return None

	hist_fid.close()

	return True

def rearrange_for_5fold_cv(lst):
	'''Puts the test fold at the end of the list.

	Messy implementation, but whatever'''

	if THIS_FOLD_OUT_OF_FIVE not in [1,2,3,4,5]:
		raise ValueError('Invalid CV fold {}'.format(THIS_FOLD_OUT_OF_FIVE))

	N = len(lst) / 5
	# print('REARRANGING FOR CV FOLD {}'.format(THIS_FOLD_OUT_OF_FIVE))
	chunks = [lst[i:i+N] for i in range(0, len(lst), N)]
	new_lst = list(itertools.chain.from_iterable(chunks[:(THIS_FOLD_OUT_OF_FIVE-1)])) + \
		list(itertools.chain.from_iterable(chunks[THIS_FOLD_OUT_OF_FIVE:])) + \
		list(itertools.chain.from_iterable(chunks[(THIS_FOLD_OUT_OF_FIVE-1):THIS_FOLD_OUT_OF_FIVE]))
	if type(lst) == type(np.array(1)):
		return np.array(new_lst)
	return new_lst
	
def get_x_data(fpath, z):
	'''
	Reads the candidate edits and returns an input tensor:
	(N_examples x N_candidate x N_maxedits x N_features)
	for each of h_lost, h_gain, bond_lost, bond_gain.

	Also needs "z", the list of true reaction strings, in order to
	generate molecules and calculate edit representations
	'''

	# Check if we've already populated these matrices once
	if os.path.isfile(fpath + '_processed'):
		with open(fpath + '_processed', 'rb') as infile:
			(x_h_lost, x_h_gain, x_bond_lost, x_bond_gain) = pickle.load(infile)	
		return (rearrange_for_5fold_cv(x_h_lost), 
			    rearrange_for_5fold_cv(x_h_gain), 
			    rearrange_for_5fold_cv(x_bond_lost), 
			    rearrange_for_5fold_cv(x_bond_gain))

	# Otherwise, we need to do it...
	with open(fpath, 'rb') as infile:
		all_candidate_edits = pickle.load(infile)
		N = len(z)
		x_h_lost = np.zeros((N, N_c, N_e, F_atom))
		x_h_gain = np.zeros((N, N_c, N_e, F_atom))
		x_bond_lost = np.zeros((N, N_c, N_e, F_bond))
		x_bond_gain = np.zeros((N, N_c, N_e, F_bond))
		for (n, candidates) in enumerate(all_candidate_edits):
			mol = Chem.MolFromSmiles(z[n].split('>')[0])
			if not mol:
				print('COULD NOT LOAD MOL: {}'.format(z[n].split('>')[0]))
			for (c, edits) in enumerate(candidates):
				if any([len(edit) > 5 for edit in edits]):
					#print('Edit counts: {}'.format([len(edit) for edit in edits]))
					#print('skipping')
					continue
				edit_h_lost_vec, edit_h_gain_vec, \
					edit_bond_lost_vec, edit_bond_gain_vec = edits_to_vectors(edits, mol)
				for (e, edit_h_lost) in enumerate(edit_h_lost_vec):
					x_h_lost[n, c, e, :] = edit_h_lost
				for (e, edit_h_gain) in enumerate(edit_h_gain_vec):
					x_h_gain[n, c, e, :] = edit_h_gain
				for (e, edit_bond_lost) in enumerate(edit_bond_lost_vec):
					x_bond_lost[n, c, e, :] = edit_bond_lost
				for (e, edit_bond_gain) in enumerate(edit_bond_gain_vec):
					x_bond_gain[n, c, e, :] = edit_bond_gain

		# Get rid of NaNs
		x_h_lost[np.isnan(x_h_lost)] = 0.0
		x_h_gain[np.isnan(x_h_gain)] = 0.0
		x_bond_lost[np.isnan(x_bond_lost)] = 0.0
		x_bond_gain[np.isnan(x_bond_gain)] = 0.0
		x_h_lost[np.isinf(x_h_lost)] = 0.0
		x_h_gain[np.isinf(x_h_gain)] = 0.0
		x_bond_lost[np.isinf(x_bond_lost)] = 0.0
		x_bond_gain[np.isinf(x_bond_gain)] = 0.0

		# Dump file so we don't have to do that again
		with open(fpath + '_processed', 'wb') as outfile:
			pickle.dump((x_h_lost, x_h_gain, x_bond_lost, x_bond_gain), outfile, pickle.HIGHEST_PROTOCOL)
		print('Converted {} to features for the first (and only) time'.format(fpath))
		return (rearrange_for_5fold_cv(x_h_lost), 
			    rearrange_for_5fold_cv(x_h_gain), 
			    rearrange_for_5fold_cv(x_bond_lost), 
			    rearrange_for_5fold_cv(x_bond_gain))

def get_xc_data(fpath):
	with open(fpath, 'rb') as infile:
		xc = pickle.load(infile)
		# Split into reagents, solvent, temp
		return (xc[:, 7:], xc[:, 1:7], xc[:, 0]) 

def get_y_data(fpath):
	# THIS IS FOR YIELD DATA, NOT PRODUCT DATA
	# ASSUME TRUE CANDIDATE ALWAYS AT FIRST INDEX
	with open(fpath, 'rb') as infile:
		y = rearrange_for_5fold_cv(pickle.load(infile))
		return np.array(y).transpose().reshape((len(y), 1)) / 100.0

def get_z_data(fpath):
	with open(fpath, 'rb') as infile:
		z = pickle.load(infile)
		return z

def pred_histogram(model, x_files, xc_files, y_files, z_files, tag = '', split_ratio = 0.8):
	'''
	Given a trained model and a list of samples, this function creates a 
	histogram of the pseudo-probabilities assigned to the true reaction 
	examples (i.e., the correct product)
	'''

	fid = open(os.path.join(FROOT, 'yields{}.dat'.format(tag)), 'w')
	fid.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
		'reaction_smiles', 'train/val', 
		'true_edit', 'yield_true_edit', 
		'predicted_edit(or no. 2)', 'yield_predicted_edit(or no. 2)',
		'rank_true_edit'
	))

	# Calculate probabilities
	
	train_preds = []; val_preds = []
	train_trues = []; val_trues = []
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
			trueprob = pred[0] # prob assigned to true outcome
			trueprobs.append(trueprob)
			rank_true_edit = 1 + len(pred) - (ss.rankdata(pred))[0]
			
			if i < int(split_ratio * preds.shape[0]):
				dataset = 'train'
				train_preds.append(trueprob)
				train_trues.append(y[i])
				if np.argmax(pred) == 0:
					train_corr += 1
			else:
				dataset = 'val'
				val_preds.append(trueprob)
				val_trues.append(y[i])
				if np.argmax(pred) == 0:
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
				edits[0], trueprob, 
				most_likely_edit, most_likely_prob,
				rank_true_edit
			))
	fid.close()
	
	train_acc = train_corr / float(len(train_preds))
	val_acc = val_corr / float(len(val_preds))

	train_preds = np.array(train_preds)
	val_preds = np.array(val_preds)
	train_trues = np.array(train_trues)
	val_trues = np.array(val_trues)

	def round3(x):
		return int(x * 1000)/1000. 

	def scatter(array1, array2, title, label):
		array1 = array1.flatten()
		array2 = array2.flatten()
		mae = np.mean(np.abs(array1 - array2))
		mse = np.mean(np.square(array1 - array2))
		try:
			plt.clf()
			plt.scatter(array1, array2, alpha = 0.5)
			plt.plot([0, 1], [0, 1], 'r--')
			plt.xlabel('True reaction yield')
			plt.ylabel('Predicted reaction yield')
			plt.title('Parity plot of yields - {} (N={})\nMAE={}, MSE={}'.format(title, len(array1), round3(mae), round3(mse)))
			plt.axis([0, 1, 0, 1])
			plt.grid(True)
			plt.savefig(os.path.join(FROOT, label), bbox_inches = 'tight')
		except:
			pass

	scatter(train_trues, train_preds, 'TRAIN', 'training_yields{}.png'.format(tag))
	scatter(val_trues, val_preds, 'VAL', 'validation_yields{}.png'.format(tag))


# def visualize_weights(model, tag):
# 	'''Given a trained model, this function visualizes the weights in each Dense layer'''
	
# 	def vec_to_png(vec, label, title):
# 		try:
# 			print(vec.shape)
# 			fig = plt.figure(figsize=(4,4))
# 			plt.pcolor(vec)#, vmin = 0, vmax = 1, cmap = plt.get_cmap('Greens'))
# 			plt.title('{} {}'.format(title, vec.shape))
# 			cbar = plt.colorbar()
# 			#plt.gca().yaxis.set_visible(False)
# 			#plt.gca().xaxis.set_visible(False)
# 			plt.xlim([0, vec.shape[1]])
# 			plt.ylim([0, vec.shape[0]])
# 			#plt.subplots_adjust(left = 0, right = 1, top = 0.4, bottom = 0)
# 			plt.savefig(os.path.join(FROOT, label) + '.png', bbox_inches = 'tight')
# 			plt.close(fig)
# 			plt.clf()
# 		except:
# 			pass

# 	for i, layer in enumerate(model.layers):
# 		try:
# 			weights = layer.layer.W.eval()
# 			vec_to_png(weights, 'weights{}(layer{})'.format(tag, i), 'Layer {}'.format(i))
# 		except Exception as e:
# 			print('Skipping layer {}'.format(i))
# 			print(e)


if __name__ == '__main__':

	FROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'output')
	
	parser = argparse.ArgumentParser()
	parser.add_argument('--nb_epoch', type = int, default = 100,
						help = 'Number of epochs to train for, default 100')
	parser.add_argument('--batch_size', type = int, default = 40,
						help = 'Batch size, default 40')
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
	parser.add_argument('--tag', type = str, default = str(int(time.time())),
						help = 'Tag for this model')
	parser.add_argument('--retrain', type = bool, default = False,
		                help = 'Retrain with loaded weights, default False')
	parser.add_argument('--test', type = bool, default = False,
						help = 'Test model only, default False')
	parser.add_argument('--l2', type = float, default = 0.01,
						help = 'l2 regularization parameter for each Dense layer, default 0.01')
	parser.add_argument('data', type = str,
		                help = 'Data folder with data files')
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
	args = parser.parse_args()

	x_files = sorted([os.path.join(args.data, dfile) \
					for dfile in os.listdir(args.data) \
					if 'candidate_edits' in dfile and '_processed' not in dfile])
	xc_files = sorted([os.path.join(args.data, dfile) \
					for dfile in os.listdir(args.data) \
					if 'context' in dfile])
	y_files = sorted([os.path.join(args.data, dfile) \
					for dfile in os.listdir(args.data) \
					if 'yield' in dfile])
	z_files = sorted([os.path.join(args.data, dfile) \
					for dfile in os.listdir(args.data) \
					if 'reaction_string' in dfile])
	print(x_files)
	print(xc_files)
	print(y_files)
	print(z_files)

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
	N_c = 100 # number of candidate edit sets
	N_e = 5 # maximum number of edits per class
	context_weight = float(args.context_weight)
	enhancement_weight = float(args.enhancement_weight)

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
			model.compile(loss = mse_of_true, 
				optimizer = Adam(lr = lr),
				)
			
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

	hist = train(model, x_files, xc_files, y_files, z_files, tag = tag, split_ratio = 0.8)
	model.save_weights(os.path.join(FROOT, 'weights{}.h5'.format(tag)), overwrite = True) 

	pred_histogram(model, x_files, xc_files, y_files, z_files, tag = tag, split_ratio = 0.8)
	#visualize_weights(model, tag)
