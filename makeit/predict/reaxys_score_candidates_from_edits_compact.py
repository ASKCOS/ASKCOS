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
from keras.optimizers import SGD, Adam, Adadelta
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

'''

This script builds, trains, tests, and saves a context-dependent reaction 
scoring model that works with a pre-generated .h5 dataset. By using a data
generator instead of reading from multiple files, data is pre-loaded for 
upcoming minibatches, so the overall training time is substantially reduced.

'''


def build(F_atom = 1, F_bond = 1, N_e = 5, N_h1 = 100, N_h2 = 50, N_h3 = 0, N_c = 100, inner_act = 'tanh', l2v = 0.01, lr = 0.0003, N_hf = 20, context_weight = 150.0, enhancement_weight = 0.1, optimizer = Adadelta()):
	'''
	Builds the feed forward model.

	N_e:  maximum number of edits of each type
	N_h1: number of hidden nodes in first layer
	N_h2: number of hidden nodes in second layer
	N_c:  maximum number of candidates for each sample (padded length)
	inner_act: activation function 
	'''

	h_lost    = Input(shape = (N_c, N_e, F_atom), name = "H_lost")
	h_gain    = Input(shape = (N_c, N_e, F_atom), name = "H_gain")
	bond_lost = Input(shape = (N_c, N_e, F_bond), name = "bond_lost")
	bond_gain = Input(shape = (N_c, N_e, F_bond), name = "bond_gain")
	reagents  = Input(shape = (256,), name = "reagent FP") # TODO: remove hard-coded length
	solvent   = Input(shape = (6,), name = "solvent descriptors c,e,s,a,b,v")
	temp      = Input(shape = (1,), name = "temperature [C]")

	h_lost_r    = Reshape((N_c * N_e, F_atom), name = "flatten H_lost")(h_lost)
	h_gain_r    = Reshape((N_c * N_e, F_atom), name = "flatten H_gain")(h_gain)
	bond_lost_r = Reshape((N_c * N_e, F_bond), name = "flatten bond_lost")(bond_lost)
	bond_gain_r = Reshape((N_c * N_e, F_bond), name = "flatten bond_gain")(bond_gain)

	h_lost_h1    = TimeDistributed(Dense(N_h1, activation = inner_act, W_regularizer = l2(l2v)), name = "embed H_lost 1")(h_lost_r)
	h_gain_h1    = TimeDistributed(Dense(N_h1, activation = inner_act, W_regularizer = l2(l2v)), name = "embed H_gain 1")(h_gain_r)
	bond_lost_h1 = TimeDistributed(Dense(N_h1, activation = inner_act, W_regularizer = l2(l2v)), name = "embed bond_lost 1")(bond_lost_r)
	bond_gain_h1 = TimeDistributed(Dense(N_h1, activation = inner_act, W_regularizer = l2(l2v)), name = "embed bond_gain 1")(bond_gain_r)
	N_h = N_h1

	if N_h2 > 0:
		h_lost_h2    = TimeDistributed(Dense(N_h2, activation = inner_act, W_regularizer = l2(l2v)), name = "embed H_lost 2")(h_lost_h1)
		h_gain_h2    = TimeDistributed(Dense(N_h2, activation = inner_act, W_regularizer = l2(l2v)), name = "embed H_gain 2")(h_gain_h1)
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
			h_lost_h    = h_lost_h2
			h_gain_h    = h_gain_h2
			bond_lost_h = bond_lost_h2
			bond_gain_h = bond_gain_h2

	else:
		h_lost_h    = h_lost_h1
		h_gain_h    = h_gain_h1
		bond_lost_h = bond_lost_h1
		bond_gain_h = bond_gain_h1

	h_lost_r2    = Reshape((N_c, N_e, N_h), name = "expand H_lost edits")(h_lost_h)
	h_gain_r2    = Reshape((N_c, N_e, N_h), name = "expand H_gain edits")(h_gain_h)
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
	temp_rpt    = RepeatVector(N_c, name = "broadcast temperature")(temp)
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
	model.compile(loss = 'categorical_crossentropy', optimizer = optimizer, 
		metrics = ['accuracy'])

	return model

def data_generator(start_at, end_at, batch_size):
	'''This function generates batches of data from the
	pickle file since all the data can't fit in memory.

	The starting and ending indices are specified explicitly so the
	same function can be used for validation data as well

	Input tensors are generated on-the-fly so there is less I/O'''

	def bond_string_to_tuple(string):
		split = string.split('-')
		return (split[0], split[1], float(split[2]))

	filePos_start_at = -1

	# Keep returning forever and ever
	with open(DATA_FPATH, 'rb') as fid:
		while True:

			# Is this the first iteration?
			if filePos_start_at == -1:

				# Remember where data starts
				legend_data = pickle.load(fid) # first doc is legend

				# Pre-load indeces
				CANDIDATE_EDITS_COMPACT = legend_data['candidate_edits_compact']
				ATOM_DESC_DICT          = legend_data['atom_desc_dict']
				T                       = legend_data['T']
				SOLVENT                 = legend_data['solvent']
				REAGENT                 = legend_data['reagent']
				YIELD                   = legend_data['yield']
				REACTION_TRUE_ONEHOT    = legend_data['reaction_true_onehot']

				for i in range(start_at): pickle.load(fid) # throw away first ___ entries
				filePos_start_at = fid.tell()
			else:
				# Jump to beginning of the desired data
				fid.seek(filePos_start_at)

			for startIndex in range(start_at, end_at, batch_size):
				endIndex = min(startIndex + batch_size, end_at)

				docs = [pickle.load(fid) for j in range(startIndex, endIndex)]

				# Initialize numpy arrays for x_h_lost, etc.
				N = endIndex - startIndex # number of samples this batch
				# Get padding ready
				x_h_lost = np.zeros((N, N_c, N_e, F_atom))
				x_h_gain = np.zeros((N, N_c, N_e, F_atom))
				x_bond_lost = np.zeros((N, N_c, N_e, F_bond))
				x_bond_gain = np.zeros((N, N_c, N_e, F_bond))
				reaction_true_onehot = np.zeros((N, N_c))

				for i, doc in enumerate(docs):

					for (c, edit_string) in enumerate(doc[CANDIDATE_EDITS_COMPACT]):
						if c >= N_c: break
						
						edit_string_split = edit_string.split(';')
						edits = [
							[atom_string for atom_string in edit_string_split[0].split(',') if atom_string],
							[atom_string for atom_string in edit_string_split[1].split(',') if atom_string],
							[bond_string_to_tuple(bond_string) for bond_string in edit_string_split[2].split(',') if bond_string],
							[bond_string_to_tuple(bond_string) for bond_string in edit_string_split[3].split(',') if bond_string],
						]

						edit_h_lost_vec, edit_h_gain_vec, \
							edit_bond_lost_vec, edit_bond_gain_vec = edits_to_vectors(edits, None, atom_desc_dict = doc[ATOM_DESC_DICT])

						for (e, edit_h_lost) in enumerate(edit_h_lost_vec):
							if e >= N_e: break
							x_h_lost[i, c, e, :] = edit_h_lost
						for (e, edit_h_gain) in enumerate(edit_h_gain_vec):
							if e >= N_e: break
							x_h_gain[i, c, e, :] = edit_h_gain
						for (e, edit_bond_lost) in enumerate(edit_bond_lost_vec):
							if e >= N_e: break
							x_bond_lost[i, c, e, :] = edit_bond_lost
						for (e, edit_bond_gain) in enumerate(edit_bond_gain_vec):
							if e >= N_e: break
							x_bond_gain[i, c, e, :] = edit_bond_gain

					# Add truncated reaction true (eventually will not truncate)
					reaction_true_onehot[i, :min(N_c, len(doc[REACTION_TRUE_ONEHOT]))] = doc[REACTION_TRUE_ONEHOT][:N_c]
				# Get rid of NaNs
				x_h_lost[np.isnan(x_h_lost)] = 0.0
				x_h_gain[np.isnan(x_h_gain)] = 0.0
				x_bond_lost[np.isnan(x_bond_lost)] = 0.0
				x_bond_gain[np.isnan(x_bond_gain)] = 0.0
				x_h_lost[np.isinf(x_h_lost)] = 0.0
				x_h_gain[np.isinf(x_h_gain)] = 0.0
				x_bond_lost[np.isinf(x_bond_lost)] = 0.0
				x_bond_gain[np.isinf(x_bond_gain)] = 0.0

				# print('Batch {} to {}'.format(startIndex, endIndex))
				# yield (x, y) as tuple, but each one is a list
				
				yield (
					[
						x_h_lost,
						x_h_gain,
						x_bond_lost,
						x_bond_gain,
						np.array([doc[REAGENT] for doc in docs]), # reagent
						np.array([doc[SOLVENT] for doc in docs]), # solvent
						np.array([doc[T] for doc in docs]), # temperature
					],
					[
						reaction_true_onehot,
					],
				)


def label_generator(start_at, end_at, batch_size):
	'''This function generates labels to match the data generated
	by data_generator'''

	filePos_start_at = -1

	# Keep returning forever and ever
	with open(LABELS_FPATH, 'rb') as fid:
		while True:
			# Is this the first iteration?
			if filePos_start_at == -1:

				# Remember where data starts
				legend_labels = pickle.load(fid) # first doc is legend
				CANDIDATE_SMILES = legend_labels['candidate_smiles']
				CANDIDATE_EDITS  = legend_labels['candidate_edits_compact']
				REACTION_TRUE    = legend_labels['reaction_true']
				RXDID            = legend_labels['rxdid']

				for i in range(start_at): pickle.load(fid) # throw away first ___ entries
				filePos_start_at = fid.tell()

			else:
				fid.seek(filePos_start_at)

			for startIndex in range(start_at, end_at, batch_size):
				endIndex = min(startIndex + batch_size, end_at)

				docs = [pickle.load(fid) for j in range(startIndex, endIndex)]
				yield {
					'candidate_smiles': [doc[CANDIDATE_SMILES] for doc in docs],
					'candidate_edits':  [doc[CANDIDATE_EDITS] for doc in docs],
					'reaction_true':    [doc[REACTION_TRUE] for doc in docs],
					'rxdid':            [doc[RXDID] for doc in docs]
				}

			filePos_start_at = -1

def get_data():
	'''Creates a dictionary defining data generators for 
	training and validation given pickled data/label files'''

	with open(DATA_FPATH, 'rb') as fid:
		legend_data = pickle.load(fid)
	with open(LABELS_FPATH, 'rb') as fid:
		legend_labels = pickle.load(fid)

	N_samples =  legend_data['N_examples']
	N_train = int(N_samples * split_ratio[0])
	N_val = int(N_samples * split_ratio[1])
	N_test = N_samples - N_train - N_val
	print('Total number of samples: {}'.format(N_samples))
	print('Training   on {}% - {}'.format(split_ratio[0]*100, N_train))
	print('Validating on {}% - {}'.format(split_ratio[1]*100, N_val))
	print('Testing    on {}% - {}'.format((1-split_ratio[1]-split_ratio[0])*100, N_test))

	return {
		'N_samples': N_samples,
		'N_train': N_train,
		#
		'train_generator': data_generator(0, N_train, batch_size),
		'train_label_generator': label_generator(0, N_train, batch_size),
		'train_nb_samples': N_train,
		#
		'val_generator': data_generator(N_train, N_train + N_val, batch_size),
		'val_label_generator': label_generator(N_train, N_train + N_val, batch_size),
		'val_nb_samples': N_val,
		#
		'test_generator': data_generator(N_train + N_val, N_samples, batch_size),
		'test_label_generator': label_generator(N_train + N_val, N_samples, batch_size),
		'test_nb_samples': N_test,
		#
		#
		'batch_size': batch_size,
	}

def train(model):
	'''Trains the Keras model'''

	data = get_data()

	# Add additional callbacks
	from keras.callbacks import ModelCheckpoint, CSVLogger, EarlyStopping
	callbacks = [
		ModelCheckpoint(WEIGHTS_FPATH, save_weights_only = True), # save every epoch
		CSVLogger(HIST_FPATH),
		EarlyStopping(patience = 5),
	]

	try:
		hist = model.fit_generator(data['train_generator'], 
			samples_per_epoch = data['train_nb_samples'],
			nb_epoch = nb_epoch, 
			validation_data = data['val_generator'],
			nb_val_samples = data['val_nb_samples'],
			#pickle_safe = True,
			callbacks = callbacks,
			verbose = 1,
		)

	except KeyboardInterrupt:
		print('Stopped training early!')

def test(model):
	'''
	Given a trained model and a list of samples, this function tests
	the model
	'''

	print('Testing model')
	data = get_data()

	fid = open(TEST_FPATH, 'w')
	fid.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
		'reaction_smiles', 'train/val', 
		'true_edit', 'prob_true_edit', 
		'predicted_edit(or no. 2)', 'prob_predicted_edit(or no. 2)',
		'rank_true_edit', 'true_smiles', 'predicted_smiles(or no. 2)',
		'RXD_id'
	))

	def test_on_set(fid, dataset, data_generator, label_generator, num_batches):
		'''Helper function that works for both training and validation sets'''
		print('Testing on {} data'.format(dataset))
		# Need to process data using generator

		our_preds = []
		true_preds = []
		corr = 0

		for batch_num in range(num_batches):
			(x, y) = data_generator.next()
			labels = label_generator.next()
			y = y[0] # only one output, which is True/False
		
			preds = model.predict_on_batch(x)

			for i in range(preds.shape[0]): 

				edits = labels['candidate_edits'][i]
				pred = preds[i, :] 
				trueprob = pred[y[i,:] != 0][0] # prob assigned to true outcome
				rank_true_edit = 1 + len(pred) - (ss.rankdata(pred))[np.argmax(y[i,:])]
				
				true_preds.append(trueprob)
				our_preds.append(pred[np.argmax(y[i,:])])
				if np.argmax(pred) == np.argmax(y[i,:]):
					corr += 1
				
				# Get most informative labels for the highest predictions
				if rank_true_edit != 1:
					# record highest probability
					most_likely_edit_i = np.argmax(pred)
					most_likely_prob   = np.max(pred)
				else:
					# record number two prediction
					most_likely_edit_i = np.argmax(pred[pred != np.max(pred)])
					most_likely_prob   = np.max(pred[pred != np.max(pred)])

				try:
					most_likely_smiles = labels['candidate_smiles'][i][most_likely_edit_i]
					most_likely_edit   = edits[most_likely_edit_i]
				except IndexError:
					most_likely_smiles = 'no_reaction'
					most_likely_edit   = 'no_reaction'

				fid.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
					labels['reaction_true'][i], dataset, 
					edits[np.argmax(y[i,:])], trueprob, 
					most_likely_edit, most_likely_prob,
					rank_true_edit, labels['reaction_true'][i].split('>')[-1], 
					most_likely_smiles, labels['rxdid'][i],
				))

		return our_preds, corr

	train_preds, train_corr = test_on_set(fid, 'train', data['train_generator'], 
		data['train_label_generator'], 
		int(np.ceil(data['train_nb_samples']/float(data['batch_size'])))
	)
	val_preds, val_corr = test_on_set(fid, 'val', data['val_generator'], 
		data['val_label_generator'], 
		int(np.ceil(data['val_nb_samples']/float(data['batch_size'])))
	)
	test_preds, test_corr = test_on_set(fid, 'test', data['test_generator'], 
		data['test_label_generator'], 
		int(np.ceil(data['test_nb_samples']/float(data['batch_size'])))
	)

	fid.close()
	
	train_acc = train_corr / float(len(train_preds))
	val_acc = val_corr / float(len(val_preds))
	test_acc = test_corr / float(len(test_preds))

	train_preds = np.array(train_preds)
	val_preds = np.array(val_preds)
	test_preds = np.array(test_preds)

	def histogram(array, title, path, acc):
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
			plt.savefig(path, bbox_inches = 'tight')
		except:
			pass

	histogram(train_preds, 'TRAIN', HISTOGRAM_FPATH % 'train', train_acc)
	histogram(val_preds, 'VAL', HISTOGRAM_FPATH % 'val', val_acc)
	histogram(test_preds, 'TEST', HISTOGRAM_FPATH % 'test', test_acc)


if __name__ == '__main__':
	
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
	parser.add_argument('--data_tag', type = str, default = 'makeit/predict/data_edits_reaxys/reaxys',
		                help = 'Data file tag, default makeit/predict/data_edits_reaxys/reaxys')
	parser.add_argument('--lr', type = float, default = 0.01, 
						help = 'Learning rate, default 0.01')
	# parser.add_argument('--dr', type = float, default = 0.5,
	# 					help = 'Dropout rate, default 0.5')
	# parser.add_argument('--fold', type = int, default = 5, 
	# 					help = 'Which fold of the 5-fold CV is this? Defaults 5')
	parser.add_argument('--visualize', type = bool, default = False,
				help = 'Whether or not to visualize weights ONLY, default False')
	parser.add_argument('--context_weight', type = float, default = 100.0,
					help = 'Weight assigned to contextual effects, default 100.0')
	parser.add_argument('--enhancement_weight', type = float, default = 0.1,
			help = 'Weight assigned to enhancement factor, default 0.1')
	parser.add_argument('--Nc', type = int, default = 500,
			help = 'Number of candidates per example, default 500')
	parser.add_argument('--optimizer', type = str, default = 'adadelta',
			help = 'Optimizer to use, default adadelta')

	args = parser.parse_args()

	mol = Chem.MolFromSmiles('[C:1][C:2]')
	(a, _, b, _) = edits_to_vectors((['1'],[],[('1','2',1.0)],[]), mol)

	F_atom = len(a[0])
	F_bond = len(b[0])
	
	nb_epoch           = int(args.nb_epoch)
	batch_size         = int(args.batch_size)
	N_h1               = int(args.Nh1)
	N_h2               = int(args.Nh2)
	N_h3               = int(args.Nh3)
	N_hf               = int(args.Nhf)
	l2v                = float(args.l2)
	lr                 = float(args.lr)
	N_c                = int(args.Nc) # number of candidate edit sets
	N_e                = 5 # maximum number of edits per class
	context_weight     = float(args.context_weight)
	enhancement_weight = float(args.enhancement_weight)
	optimizer          = args.optimizer

	# THIS_FOLD_OUT_OF_FIVE = int(args.fold)
	tag = args.tag

	split_ratio = (0.8, 0.1) # 80% training, 10% validation, balance testing

	if optimizer == 'sgd':
		opt  = SGD(lr = lr, decay = 1e-4, momentum = 0.9)
	elif optimizer == 'adam':
		opt = Adam(lr = lr)
	elif optimizer == 'adadelta':
		opt = Adadelta()
		print('Because Adadelta was selected, ignoring lr setting')
	else:
		raise ValueError('Unrecognized optimizer')

	# Labels
	FROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'output')
	FROOT = os.path.join(FROOT, tag)
	if not os.path.isdir(FROOT):
		os.mkdir(FROOT)
	MODEL_FPATH = os.path.join(FROOT, '{} model.json'.format(tag))
	WEIGHTS_FPATH = os.path.join(FROOT, '{} weights.h5'.format(tag))
	HIST_FPATH = os.path.join(FROOT, '{} hist.csv'.format(tag))
	TEST_FPATH = os.path.join(FROOT, '{} probs.dat'.format(tag))
	HISTOGRAM_FPATH = os.path.join(FROOT, '{} histogram %s.png'.format(tag))

	DATA_FPATH = '{}_data.pickle'.format(args.data_tag)
	LABELS_FPATH = '{}_labels.pickle'.format(args.data_tag)

	if bool(args.retrain):
		print('Reloading from file')
		rebuild = raw_input('Do you want to rebuild from scratch instead of loading from file? [n/y] ')
		if rebuild == 'y':
			model = build(F_atom = F_atom, F_bond = F_bond, N_e = N_e, N_c = N_c, N_h1 = N_h1, 
				N_h2 = N_h2, N_h3 = N_h3, N_hf = N_hf, l2v = l2v, lr = lr, optimizer = opt,
				context_weight = context_weight, enhancement_weight = enhancement_weight)
		else:
			model = model_from_json(open(MODEL_FPATH).read())
			model.compile(loss = 'categorical_crossentropy', 
				optimizer = opt,
				metrics = ['accuracy'])
		model.load_weights(WEIGHTS_FPATH)
	else:
		model = build(F_atom = F_atom, F_bond = F_bond, N_e = N_e, N_c = N_c, N_h1 = N_h1, N_h2 = N_h2, N_h3 = N_h3, N_hf = N_hf, l2v = l2v, lr = lr, context_weight = context_weight,
			enhancement_weight = enhancement_weight, optimizer = opt)
		try:
			with open(MODEL_FPATH, 'w') as outfile:
				outfile.write(model.to_json())
		except:
			print('could not write model to json')

	if bool(args.test):
		test(model)
		quit(1)

	train(model)
	model.save_weights(WEIGHTS_FPATH, overwrite = True) 
	test(model)