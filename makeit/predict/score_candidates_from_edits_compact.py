# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import numpy as np
import os
import sys
import argparse
import h5py # needed for save_weights, fails otherwise
from keras import backend as K 
from keras.models import Sequential, Model, model_from_json
from keras.layers import Dense, Activation, Input
from keras.layers.core import Flatten, Permute, Reshape, Dropout, Lambda
from keras.layers.wrappers import TimeDistributed
from keras.engine.topology import Merge, merge
from keras.optimizers import *
from keras.layers.convolutional import Convolution1D, Convolution2D
from keras.regularizers import l2
from keras.utils.np_utils import to_categorical
from makeit.embedding.descriptors import edits_to_vectors, oneHotVector # for testing
import rdkit.Chem as Chem
import theano.tensor as T
import cPickle as pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt    # for visualization
import scipy.stats as ss
import itertools
import time
from makeit.utils.threadsafe import threadsafe_generator


def build(F_atom = 1, F_bond = 1, N_h1 = 100, N_h2 = 50, N_h3 = 0, inner_act = 'tanh', l2v = 0.0, lr = 0.0003, N_hf = 20, optimizer = Adadelta(), extra_outputs = False, absolute_score = False):
	'''
	Builds the feed forward model.

	N_e:  maximum number of edits of each type
	N_h1: number of hidden nodes in first layer
	N_h2: number of hidden nodes in second layer
	inner_act: activation function 
	'''
	if BASELINE_MODEL:
		FPs = Input(shape = (None, 1024), name = "FPs")
		features = TimeDistributed(Dense(N_hf, activation = inner_act), name = "FPs to features")(FPs)
		unscaled_score = TimeDistributed(Dense(1, activation = 'linear'), name = "features to score")(features)
		dynamic_flattener       = lambda x: T.reshape(x, (x.shape[0], x.shape[1]), ndim  = x.ndim-1)
		dynamic_flattener_shape = lambda x: (None, x[1])
		unscaled_score_flat    = Lambda(dynamic_flattener, output_shape = dynamic_flattener_shape, name = "flatten")(unscaled_score)
		if absolute_score:
			score = unscaled_score_flat 
		else:
			score = Activation('softmax', name = "scores to probs")(unscaled_score_flat)
		model = Model(input = [FPs], 
			output = [score])
		model.summary()
		model.compile(loss = 'categorical_crossentropy', optimizer = optimizer, 
				metrics = ['accuracy'])
		return model

	h_lost    = Input(shape = (None, None, F_atom), name = "H_lost")
	h_gain    = Input(shape = (None, None, F_atom), name = "H_gain")
	bond_lost = Input(shape = (None, None, F_bond), name = "bond_lost")
	bond_gain = Input(shape = (None, None, F_bond), name = "bond_gain")

	# Combine along first three dimensions
	dynamic_reshaper       = lambda x: T.reshape(x, (x.shape[0] * x.shape[1] * x.shape[2], x.shape[3]), ndim  = x.ndim-2)
	dynamic_reshaper_shape = lambda x: (None,) + x[3:]

	h_lost_r    = Lambda(dynamic_reshaper, output_shape = dynamic_reshaper_shape, name = "flatten_H_lost")(h_lost)
	h_gain_r    = Lambda(dynamic_reshaper, output_shape = dynamic_reshaper_shape, name = "flatten_H_gain")(h_gain)
	bond_lost_r = Lambda(dynamic_reshaper, output_shape = dynamic_reshaper_shape, name = "flatten_bond_lost")(bond_lost)
	bond_gain_r = Lambda(dynamic_reshaper, output_shape = dynamic_reshaper_shape, name = "flatten_bond_gain")(bond_gain)

	h_lost_h1    = Dense(N_h1, activation = inner_act, W_regularizer = l2(l2v), name = "embed H_lost 1")(h_lost_r)
	h_gain_h1    = Dense(N_h1, activation = inner_act, W_regularizer = l2(l2v), name = "embed H_gain 1")(h_gain_r)
	bond_lost_h1 = Dense(N_h1, activation = inner_act, W_regularizer = l2(l2v), name = "embed bond_lost 1")(bond_lost_r)
	bond_gain_h1 = Dense(N_h1, activation = inner_act, W_regularizer = l2(l2v), name = "embed bond_gain 1")(bond_gain_r)

	N_h = N_h1

	if N_h2 > 0:
		h_lost_h2    = Dense(N_h2, activation = inner_act, W_regularizer = l2(l2v), name = "embed H_lost 2")(h_lost_h1)
		h_gain_h2    = Dense(N_h2, activation = inner_act, W_regularizer = l2(l2v), name = "embed H_gain 2")(h_gain_h1)
		bond_lost_h2 = Dense(N_h2, activation = inner_act, W_regularizer = l2(l2v), name = "embed bond_lost 2")(bond_lost_h1)
		bond_gain_h2 = Dense(N_h2, activation = inner_act, W_regularizer = l2(l2v), name = "embed bond_gain 2")(bond_gain_h1)
		N_h = N_h2

		if N_h3 > 0:
			h_lost_h    = Dense(N_h3, activation = inner_act, W_regularizer = l2(l2v), name = "embed H_lost 3")(h_lost_h2)
			h_gain_h    = Dense(N_h3, activation = inner_act, W_regularizer = l2(l2v), name = "embed H_gain 3")(h_gain_h2)
			bond_lost_h = Dense(N_h3, activation = inner_act, W_regularizer = l2(l2v), name = "embed bond_lost 3")(bond_lost_h2)
			bond_gain_h = Dense(N_h3, activation = inner_act, W_regularizer = l2(l2v), name = "embed bond_gain 3")(bond_gain_h2)
			N_h         = N_h3

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

	# Re-expand (using tricky Merge layer, where x[0] is actual data and x[1] is only used for shape)
	dynamic_unreshaper = lambda x: T.reshape(x[0], (x[1].shape[0], x[1].shape[1], x[1].shape[2], x[0].shape[1]), ndim  = x[0].ndim+2)
	dynamic_unreshaper_shape = lambda x: x[1][:3] + x[0][1:]

	h_lost_r2    = merge([h_lost_h, h_lost], mode = dynamic_unreshaper, output_shape = dynamic_unreshaper_shape, name = "expand H_lost edits")
	h_gain_r2    = merge([h_gain_h, h_gain], mode = dynamic_unreshaper, output_shape = dynamic_unreshaper_shape, name = "expand H_gain edits")
	bond_lost_r2 = merge([bond_lost_h, bond_lost], mode = dynamic_unreshaper, output_shape = dynamic_unreshaper_shape, name = "expand bond_lost edits")
	bond_gain_r2 = merge([bond_gain_h, bond_gain], mode = dynamic_unreshaper, output_shape = dynamic_unreshaper_shape, name = "expand bond_gain edits")

	# Add edits within a single candidate
	sum_along_axis2       = lambda x: K.sum(x, axis = 2)
	sum_along_axis2_shape = lambda x: x[:2] + x[3:]
	h_lost_sum    = Lambda(sum_along_axis2, output_shape = sum_along_axis2_shape, name = "sum H_lost")(h_lost_r2)
	h_gain_sum    = Lambda(sum_along_axis2, output_shape = sum_along_axis2_shape, name = "sum H_gain")(h_gain_r2)
	bond_lost_sum = Lambda(sum_along_axis2, output_shape = sum_along_axis2_shape, name = "sum bond_lost")(bond_lost_r2)
	bond_gain_sum = Lambda(sum_along_axis2, output_shape = sum_along_axis2_shape, name = "sum bond_gain")(bond_gain_r2)

	# Sum across edits in their intermediate representation
	net_sum = merge([h_lost_sum, h_gain_sum, bond_lost_sum, bond_gain_sum], mode = 'sum', name = "sum_across edits")

	feature_to_feature = Dense(N_hf, activation = inner_act, W_regularizer = l2(l2v))
	net_sum_h = TimeDistributed(feature_to_feature, name = "reaction embedding post-sum")(net_sum)

	unscaled_score = TimeDistributed(Dense(1, activation = 'linear'), name = "feature to score")(net_sum_h)

	# Flatten
	# Combine along first three dimensions
	dynamic_flattener       = lambda x: T.reshape(x, (x.shape[0], x.shape[1]), ndim  = x.ndim-1)
	dynamic_flattener_shape = lambda x: (None, x[1])

	unscaled_score_flat    = Lambda(dynamic_flattener, output_shape = dynamic_flattener_shape, name = "flatten")(unscaled_score)

	if absolute_score:
		score = unscaled_score_flat 
	else:
		score = Activation('softmax', name = "scores to probs")(unscaled_score_flat)
	
	#score = unscaled_score_r

	if extra_outputs:
		model = Model(input = [h_lost, h_gain, bond_lost, bond_gain, reagents, solvent, temp], 
			output = [h_lost_sum, h_gain_sum, bond_lost_sum, bond_gain_sum, net_sum, net_sum_h, params, unscaled_score, score])
		return model

	model = Model(input = [h_lost, h_gain, bond_lost, bond_gain], 
		output = [score])

	model.summary()

	# Now compile
	model.compile(loss = 'categorical_crossentropy', optimizer = optimizer, 
			metrics = ['accuracy'])
	
	return model

@threadsafe_generator
def data_generator(start_at, end_at, batch_size, max_N_c = None, shuffle = False, allowable_batchNums = set()):
	'''This function generates batches of data from the
	pickle file since all the data can't fit in memory.

	The starting and ending indices are specified explicitly so the
	same function can be used for validation data as well

	Input tensors are generated on-the-fly so there is less I/O

	max_N_c is the maximum number of candidates to consider. This should ONLY be used
	for training, not for validation or testing.

	"mybatchnums" is a new list that contains the batch indices (across the whole dataset) that belong
	to this particular generator. This allows for CV splitting *outside* of this function.'''

	def bond_string_to_tuple(string):
		split = string.split('-')
		return (split[0], split[1], float(split[2]))

	fileInfo  = [() for j in range(start_at, end_at, batch_size)] # (filePos, startIndex, endIndex)
	batchDims = [() for j in range(start_at, end_at, batch_size)] # dimensions of each batch
	batchNums = np.array([i for (i, j) in enumerate(range(start_at, end_at, batch_size))]) # list to shuffle later

	# Keep returning forever and ever
	with open(DATA_FPATH, 'rb') as fid:

		# Do a first pass through the data
		legend_data = pickle.load(fid) # first doc is legend

		# Pre-load indeces
		CANDIDATE_EDITS_COMPACT = legend_data['candidate_edits_compact']
		ATOM_DESC_DICT          = legend_data['atom_desc_dict']
		REACTION_TRUE_ONEHOT    = legend_data['reaction_true_onehot']

		for i in range(start_at): pickle.load(fid) # throw away first ___ entries

		for k, startIndex in enumerate(range(start_at, end_at, batch_size)):
			endIndex = min(startIndex + batch_size, end_at)

			# Remember this starting position
			fileInfo[k] = (fid.tell(), startIndex, endIndex)

			N = endIndex - startIndex # number of samples this batch
			# print('Serving up examples {} through {}'.format(startIndex, endIndex))

			docs = [pickle.load(fid) for j in range(startIndex, endIndex)]

			# FNeed to figure out size of padded batch
			N_c = max([len(doc[REACTION_TRUE_ONEHOT]) for doc in docs])
			if type(max_N_c) != type(None): # allow truncation during training
				N_c = min(N_c, max_N_c)
			N_e1 = 1; N_e2 = 1; N_e3 = 1; N_e4 = 1;
			for i, doc in enumerate(docs):
				for (c, edit_string) in enumerate(doc[CANDIDATE_EDITS_COMPACT]):
					if c >= N_c: break
					edit_string_split = edit_string.split(';')
					N_e1 = max(N_e1, edit_string_split[0].count(',') + 1)
					N_e2 = max(N_e2, edit_string_split[1].count(',') + 1)
					N_e3 = max(N_e3, edit_string_split[2].count(',') + 1)
					N_e4 = max(N_e4, edit_string_split[3].count(',') + 1)

			# Remember sizes of x_h_lost, x_h_gain, x_bond_lost, x_bond_gain, reaction_true_onehot
			batchDim = (N, N_c, N_e1, N_e2, N_e3, N_e4)

			# print('The padded sizes of this batch will be: N, N_c, N_e1, N_e2, N_e3, N_e4')
			# print(batchDim)
			batchDims[k] = batchDim

		while True:

			if shuffle: np.random.shuffle(batchNums)

			for batchNum in batchNums:
				if batchNum not in allowable_batchNums: continue

				(filePos, startIndex, endIndex) = fileInfo[batchNum]
				(N, N_c, N_e1, N_e2, N_e3, N_e4) = batchDims[batchNum]
				fid.seek(filePos)

				N = endIndex - startIndex # number of samples this batch
				# print('Serving up examples {} through {}'.format(startIndex, endIndex))

				docs = [pickle.load(fid) for j in range(startIndex, endIndex)]

				if BASELINE_MODEL:
					x = np.zeros((N, N_c, 1024), dtype = np.float32) 

				# Initialize numpy arrays for x_h_lost, etc.
				x_h_lost = np.zeros((N, N_c, N_e1, F_atom), dtype=np.float32)
				x_h_gain = np.zeros((N, N_c, N_e2, F_atom), dtype=np.float32)
				x_bond_lost = np.zeros((N, N_c, N_e3, F_bond), dtype=np.float32)
				x_bond_gain = np.zeros((N, N_c, N_e4, F_bond), dtype=np.float32)
				reaction_true_onehot = np.zeros((N, N_c), dtype=np.float32)

				for i, doc in enumerate(docs):

					for (c, edit_string) in enumerate(doc[CANDIDATE_EDITS_COMPACT]):
						if c >= N_c: 
							break

						if BASELINE_MODEL:
							x[i, c, :] = doc[legend_data['prod_FPs']][c]
						
						edit_string_split = edit_string.split(';')
						edits = [
							[atom_string for atom_string in edit_string_split[0].split(',') if atom_string],
							[atom_string for atom_string in edit_string_split[1].split(',') if atom_string],
							[bond_string_to_tuple(bond_string) for bond_string in edit_string_split[2].split(',') if bond_string],
							[bond_string_to_tuple(bond_string) for bond_string in edit_string_split[3].split(',') if bond_string],
						]

						try:
							edit_h_lost_vec, edit_h_gain_vec, \
								edit_bond_lost_vec, edit_bond_gain_vec = edits_to_vectors(edits, None, atom_desc_dict = doc[ATOM_DESC_DICT], ORIGINAL_VERSION = True)
						except KeyError as e: # sometimes molAtomMapNumber not found if hydrogens were explicit
							continue

						for (e, edit_h_lost) in enumerate(edit_h_lost_vec):
							if e >= N_e1: raise ValueError('N_e1 not large enough!')
							x_h_lost[i, c, e, :] = edit_h_lost
						for (e, edit_h_gain) in enumerate(edit_h_gain_vec):
							if e >= N_e2: raise ValueError('N_e2 not large enough!')
							x_h_gain[i, c, e, :] = edit_h_gain
						for (e, edit_bond_lost) in enumerate(edit_bond_lost_vec):
							if e >= N_e3: raise ValueError('N_e3 not large enough!')
							x_bond_lost[i, c, e, :] = edit_bond_lost
						for (e, edit_bond_gain) in enumerate(edit_bond_gain_vec):
							if e >= N_e4: raise ValueRrror('N_e4 not large enough!')
							x_bond_gain[i, c, e, :] = edit_bond_gain

					# Add truncated reaction true (eventually will not truncate)
					if type(max_N_c) == type(None):
						reaction_true_onehot[i, :len(doc[REACTION_TRUE_ONEHOT])] = doc[REACTION_TRUE_ONEHOT]
					else:
						reaction_true_onehot[i, :min(len(doc[REACTION_TRUE_ONEHOT]), max_N_c)] = doc[REACTION_TRUE_ONEHOT][:max_N_c]

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

				y = reaction_true_onehot

				if BASELINE_MODEL:
					yield ([x], [y])

				else:

					yield (
						[
							x_h_lost,
							x_h_gain,
							x_bond_lost,
							x_bond_gain,
						],
						[
							y,
						],
					)

@threadsafe_generator
def label_generator(start_at, end_at, batch_size, allowable_batchNums = set()):
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

			for batchNum, startIndex in enumerate(range(start_at, end_at, batch_size)):
				if batchNum not in allowable_batchNums: continue

				endIndex = min(startIndex + batch_size, end_at)

				docs = [pickle.load(fid) for j in range(startIndex, endIndex)]
				yield {
					'candidate_smiles': [doc[CANDIDATE_SMILES] for doc in docs],
					'candidate_edits':  [doc[CANDIDATE_EDITS] for doc in docs],
					'reaction_true':    [doc[REACTION_TRUE] for doc in docs],
					'rxdid':            [doc[RXDID] for doc in docs]
				}

			filePos_start_at = -1

def get_data(max_N_c = None, shuffle = False):
	'''Creates a dictionary defining data generators for 
	training and validation given pickled data/label files

	max_N_c and shuffle only refers to training data'''

	with open(DATA_FPATH, 'rb') as fid:
		legend_data = pickle.load(fid)
	with open(LABELS_FPATH, 'rb') as fid:
		legend_labels = pickle.load(fid)

	N_samples =  legend_data['N_examples']

	# New approach: each generator takes whole set, but use allowable_batchNums to filter down
	from sklearn.cross_validation import KFold
	pseudoRandomCV = KFold(len(range(0, N_samples, batch_size)), n_folds = 5, shuffle = True, random_state = 0)
	
	# Stupid solution since KFold can't be indexed
	for i, (trainval_batches, test_batches) in enumerate(pseudoRandomCV):
		if i == (FOLD_NUM - 1): break

	print(trainval_batches)
	np.random.seed(0)
	np.random.shuffle(trainval_batches)
	train_batches = trainval_batches[:int(len(trainval_batches) * split_ratio[0] / (split_ratio[0] + split_ratio[1]))]
	val_batches   = trainval_batches[int(len(trainval_batches) * split_ratio[0] / (split_ratio[0] + split_ratio[1])):]


	print('Train batches: {}'.format(train_batches))
	print('Val batches: {}'.format(val_batches))
	print('Test batches: {}'.format(test_batches))
	
	N_train = len(train_batches) * batch_size
	N_val = len(val_batches) * batch_size
	N_test = len(test_batches) * batch_size
	print('Total number of samples: {}'.format(N_samples))
	print('Training   on {}% - {}'.format(split_ratio[0]*100, N_train))
	print('Validating on {}% - {}'.format(split_ratio[1]*100, N_val))
	print('Testing    on {}% - {}'.format((1-split_ratio[1]-split_ratio[0])*100, N_test))

	# return {
	# 	'N_samples': N_samples,
	# 	'N_train': N_train,
	# 	#
	# 	'train_generator': data_generator(0, N_train, batch_size, max_N_c = max_N_c, shuffle = shuffle),
	# 	'train_label_generator': label_generator(0, N_train, batch_size),
	# 	'train_nb_samples': N_train,
	# 	#
	# 	'val_generator': data_generator(N_train, N_train + N_val, batch_size),
	# 	'val_label_generator': label_generator(N_train, N_train + N_val, batch_size),
	# 	'val_nb_samples': N_val,
	# 	#
	# 	'test_generator': data_generator(N_train + N_val, N_samples, batch_size),
	# 	'test_label_generator': label_generator(N_train + N_val, N_samples, batch_size),
	# 	'test_nb_samples': N_test,
	# 	#
	# 	#
	# 	'batch_size': batch_size,
	# }



	return {
		'N_samples': N_samples,
		'N_train': N_train,
		#
		'train_generator': data_generator(0, N_samples, batch_size, max_N_c = max_N_c, shuffle = shuffle, allowable_batchNums = train_batches),
		'train_label_generator': label_generator(0, N_samples, batch_size, allowable_batchNums = train_batches),
		'train_nb_samples': N_train,
		#
		'val_generator': data_generator(0, N_samples, batch_size, allowable_batchNums = val_batches),
		'val_label_generator': label_generator(0, N_samples, batch_size, allowable_batchNums = val_batches),
		'val_nb_samples': N_val,
		#
		'test_generator': data_generator(0, N_samples, batch_size, allowable_batchNums = test_batches),
		'test_label_generator': label_generator(0, N_samples, batch_size, allowable_batchNums = test_batches),
		'test_nb_samples': N_test,
		#
		#
		'batch_size': batch_size,
	}

def train(model, data):
	'''Trains the Keras model'''

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

def test(model, data):
	'''
	Given a trained model and a list of samples, this function tests
	the model
	'''

	print('Testing model')

	fid = open(TEST_FPATH, 'w')
	fid.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
		'reaction_smiles', 'train/val', 
		'true_edit', 'prob_true_edit', 
		'predicted_edit(or no. 2)', 'prob_predicted_edit(or no. 2)',
		'rank_true_edit', 'true_smiles', 'predicted_smiles(or no. 2)',
		'RXD_id','yield_true',
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
			y = y[0] # only one output, which is True/False or yield
		
			# TODO: pre-fetch data in queue
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
				trueyield = 0.0

				try:
					most_likely_smiles = labels['candidate_smiles'][i][most_likely_edit_i]
					most_likely_edit   = edits[most_likely_edit_i]
				except IndexError:
					most_likely_smiles = 'no_reaction'
					most_likely_edit   = 'no_reaction'


				fid.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
					labels['reaction_true'][i], dataset, 
					edits[np.argmax(y[i,:])], trueprob, 
					most_likely_edit, most_likely_prob,
					rank_true_edit, labels['reaction_true'][i].split('>')[-1], 
					most_likely_smiles, labels['rxdid'][i], trueyield
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

	histogram(train_preds, 'TRAIN', HISTOGRAM_FPATH('train'), train_acc)
	histogram(val_preds, 'VAL', HISTOGRAM_FPATH('val'), val_acc)
	histogram(test_preds, 'TEST', HISTOGRAM_FPATH('test'), test_acc)


if __name__ == '__main__':

	np.random.seed(0)
	
	parser = argparse.ArgumentParser()
	parser.add_argument('--nb_epoch', type = int, default = 200,
						help = 'Max number of epochs to train for, default 200')
	parser.add_argument('--batch_size', type = int, default = 20,
						help = 'Batch size, default 200')
	parser.add_argument('--Nh1', type = int, default = 200,
						help = 'Number of hidden nodes in first layer, default 100')
	parser.add_argument('--Nh2', type = int, default = 100,
						help = 'Number of hidden nodes in second layer, default 50')
	parser.add_argument('--Nh3', type = int, default = 50,
						help = 'Number of hidden nodes in third layer, ' + 
								'immediately before summing, default 50')
	parser.add_argument('--Nhf', type = int, default = 50,
						help = 'Number of hidden nodes in layer between summing ' +
								'and final score, default 50')
	parser.add_argument('--tag', type = str, default = str(time.time()),
						help = 'Tag for this model')
	parser.add_argument('--retrain', type = bool, default = False,
		                help = 'Retrain with loaded weights, default False')
	parser.add_argument('--test', type = bool, default = False,
						help = 'Test model only, default False')
	parser.add_argument('--l2', type = float, default = 0.0,
						help = 'l2 regularization parameter for each Dense layer, default 0.0')
	parser.add_argument('--data_tag', type = str, default = 'makeit/predict/lowe_data_edits/lowe',
		                help = 'Data file tag, default makeit/predict/lowe_data_edits/lowe')
	parser.add_argument('--lr', type = float, default = 0.01, 
						help = 'Learning rate, default 0.01')
	parser.add_argument('--Nc', type = int, default = 1000,
			help = 'Number of candidates to truncate to during training, default 1000')
	parser.add_argument('--optimizer', type = str, default = 'adadelta',
			help = 'Optimizer to use, default adadelta')
	parser.add_argument('--inner_act', type = str, default = 'tanh',
			help = 'Inner activation function, default "tanh" ')
	parser.add_argument('--fold', type = int, default = 1, 
						help = 'Which fold of the 5-fold CV is this? Defaults 1')
	parser.add_argument('--baseline', type = int, default = 0 ,
					help = 'Baseline fingerprint model? Default 0')

	args = parser.parse_args()

	
	nb_epoch           = int(args.nb_epoch)
	batch_size         = int(args.batch_size)
	N_h1               = int(args.Nh1)
	N_h2               = int(args.Nh2)
	N_h3               = int(args.Nh3)
	N_hf               = int(args.Nhf)
	l2v                = float(args.l2)
	lr                 = float(args.lr)
	max_N_c            = int(args.Nc) # number of candidate edit sets
	optimizer          = args.optimizer
	inner_act          = args.inner_act
	FOLD_NUM           = int(args.fold)
	BASELINE_MODEL     = bool(int(args.baseline))

	# THIS_FOLD_OUT_OF_FIVE = int(args.fold)
	tag = args.tag + str(args.fold)

	split_ratio = (0.7, 0.1) # 80% training, 10% validation, balance testing

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
	FROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
	FROOT = os.path.join(FROOT, tag)
	if not os.path.isdir(FROOT):
		print(FROOT)
		os.mkdir(FROOT)
	MODEL_FPATH = os.path.join(FROOT, 'model.json')
	WEIGHTS_FPATH = os.path.join(FROOT, 'weights.h5')
	HIST_FPATH = os.path.join(FROOT, 'hist.csv')
	TEST_FPATH = os.path.join(FROOT, 'probs.dat')
	HISTOGRAM_FPATH = lambda x: os.path.join(FROOT, 'histogram {}.png'.format(x))
	ARGS_FPATH = os.path.join(FROOT, 'args.json')

	with open(ARGS_FPATH, 'w') as fid:
		import json
		json.dump(args.__dict__, fid)

	DATA_FPATH = '{}_data.pickle'.format(args.data_tag)
	LABELS_FPATH = '{}_labels.pickle'.format(args.data_tag)

	this_dir = os.getcwd()
	mol = Chem.MolFromSmiles('[CH3:1][CH3:2]')
	(a, _, b, _) = edits_to_vectors((['1'],[],[('1','2',1.0)],[]), mol, ORIGINAL_VERSION = True)
	os.chdir(this_dir)

	F_atom = len(a[0])
	F_bond = len(b[0])

	if bool(args.retrain):
		print('Reloading from file')
		rebuild = raw_input('Do you want to rebuild from scratch instead of loading from file? [n/y] ')
		if rebuild == 'y':
			model = build(F_atom = F_atom, F_bond = F_bond, N_h1 = N_h1, 
				N_h2 = N_h2, N_h3 = N_h3, N_hf = N_hf, l2v = l2v, lr = lr, optimizer = opt, inner_act = inner_act,
				extra_outputs = bool(args.visualize))
		else:
			model = model_from_json(open(MODEL_FPATH).read())
			model.compile(loss = 'categorical_crossentropy', 
			optimizer = opt,
			metrics = ['accuracy'])
		model.load_weights(WEIGHTS_FPATH)
	else:
		model = build(F_atom = F_atom, F_bond = F_bond, N_h1 = N_h1, N_h2 = N_h2, N_h3 = N_h3, N_hf = N_hf, l2v = l2v, lr = lr, inner_act = inner_act, optimizer = opt)
		try:
			with open(MODEL_FPATH, 'w') as outfile:
				outfile.write(model.to_json())
		except:
			print('could not write model to json')

	if bool(args.test):
		data = get_data(max_N_c = max_N_c, shuffle = False)
		test(model, data)
		quit(1)
	data = get_data(max_N_c = max_N_c, shuffle = True)
	train(model, data)
	model.save_weights(WEIGHTS_FPATH, overwrite = True) 
	data = get_data(max_N_c = max_N_c, shuffle = False)

	test(model, data)