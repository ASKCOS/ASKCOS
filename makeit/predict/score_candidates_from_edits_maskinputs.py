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

def get_y_data(fpath):
	with open(fpath, 'rb') as infile:
		y = rearrange_for_5fold_cv(pickle.load(infile))
		y = np.array([y_i.index(True) for y_i in y])
		return to_categorical(y, nb_classes = N_c)

def get_z_data(fpath):
	# DONT REARRANGE FOR CV UNTIL AFTER GETTING X DATA
	with open(fpath, 'rb') as infile:
		z = pickle.load(infile)
		return z

def get_accuracy(model, x_files, y_files, z_files, tag = '', split_ratio = 0.8, cripple = None):
	'''
	Given a trained model and a list of samples, this function creates a 
	histogram of the pseudo-probabilities assigned to the true reaction 
	examples (i.e., the correct product)
	'''

	# Calculate probabilities
	
	train_preds = []; val_preds = []
	train_corr = 0; val_corr = 0;
	dataset = ''; trueprob = 0; trueprobs = []

	for fnum in range(len(x_files)):
		print('Testing file number {}'.format(fnum))
		# Data must be pre-padded

		z = get_z_data(z_files[fnum])
		(x_h_lost, x_h_gain, x_bond_lost, x_bond_gain) = list(get_x_data(x_files[fnum], z))
		y = get_y_data(y_files[fnum])
		z = rearrange_for_5fold_cv(z)

		if cripple is not None:
			for i in cripple:
				offset = x_bond_lost.shape[-1] - x_h_lost.shape[-1]
				# Set atom attributes to zero
				x_h_lost[:, :, :, i] = 0 # H_lost
				x_h_gain[:, :, :, i] = 0 # H_gain
				# Set atom attributes within bond attributes to zero
				x_bond_lost[:, :, :, i] = 0 # Bond_lost
				x_bond_lost[:, :, :, i + offset] = 0 # Bond_lost
				x_bond_gain[:, :, :, i] = 0 # Bond_gain
				x_bond_gain[:, :, :, i + offset] = 0 # Bond_gain
				print('Set index {} and {} to zero'.format(i, i + offset))
		
		# Getting unprocessed x data
		with open(x_files[fnum], 'rb') as infile:
			all_edits = pickle.load(infile)
                all_edits = rearrange_for_5fold_cv(all_edits)

		preds = model.predict([x_h_lost, x_h_gain, x_bond_lost, x_bond_gain], batch_size = 20)
		print(preds)

		for i in range(preds.shape[0]): 
			rank_true_edit = 0
			pred = preds[i, :] # Iterate through each sample
			trueprob = pred[y[i,:] != 0][0] # prob assigned to true outcome
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
	
	train_acc = train_corr / float(len(train_preds))
	val_acc = val_corr / float(len(val_preds))

	return (train_acc, val_acc)
	

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
	parser.add_argument('--Nhf', type = int, default = 0,
						help = 'Number of hidden nodes in layer between summing ' +
								'and final score, default 0')
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
	parser.add_argument('--Nc', type = int, default = 500,
						help = 'Number of candidates per example, default 500')
	# parser.add_argument('--dr', type = float, default = 0.5,
	# 					help = 'Dropout rate, default 0.5')
	parser.add_argument('--visualize', type = bool, default = False,
		                help = 'Whether or not to visualize weights ONLY, default False')
	parser.add_argument('--fold', type = int, default = 5, 
						help = 'Which fold of the 5-fold CV is this? Defaults 5')
	parser.add_argument('--optimizer', type = str, default = 'SGD',
						help = 'Optimizer? Adam or SGD right now, default SGD')
	args = parser.parse_args()

	x_files = sorted([os.path.join(args.data, dfile) \
					for dfile in os.listdir(args.data) \
					if 'candidate_edits' in dfile and '_processed' not in dfile])
	y_files = sorted([os.path.join(args.data, dfile) \
					for dfile in os.listdir(args.data) \
					if 'candidate_bools' in dfile])
	z_files = sorted([os.path.join(args.data, dfile) \
					for dfile in os.listdir(args.data) \
					if 'reaction_string' in dfile])
	print(x_files)
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
	N_c =  500 # number of candidate edit sets
	N_e = 5 # maximum number of edits per class

	THIS_FOLD_OUT_OF_FIVE = int(args.fold)
	tag = args.tag + ' fold{}'.format(args.fold)

	print('Reloading from file')
	model = model_from_json(open(os.path.join(FROOT, 'model{}.json'.format(tag))).read())
	model.compile(loss = 'categorical_crossentropy',
     	optimizer = SGD(lr = lr),
		metrics = ['accuracy']
	)
	model.load_weights(os.path.join(FROOT, 'weights{}.h5'.format(tag)), by_name = True)
	print('Loaded weights from file')
	
	fid = open(os.path.join(FROOT, 'input_masking_{}.csv'.format(tag)), 'w')
	fid.write('crippled_index\ttrain_acc\tval_acc\n')

	cripple_list = [
		[0],
		[1],
		[2],
		[3],
		[4],
		[5],
		[6],
		range(7, 18),
		range(18, 24),
		range(24, 29),
		[29],
		[30],
		[31]
	]
	for cripple in cripple_list:
		print('CRIPPLING {}'.format(cripple))
		# Cripple index i
		(train_acc, val_acc) = get_accuracy(model, x_files, y_files, z_files, tag = tag, split_ratio = 0.8, cripple = cripple)
		fid.write('{}\t{}\t{}\n'.format(cripple, train_acc, val_acc))
