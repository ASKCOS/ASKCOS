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
from makeit.predict.preprocess_candidates import *
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
		x = list(get_x_data(x_files[fnum], z))
		y = get_y_data(y_files[fnum])
		z = rearrange_for_5fold_cv(z)

		if cripple:
			offset = x[2].shape[-1] - x[0].shape[-1]
			# Set atom attributes to zero
			x[0][:, :, :, cripple] = 0 # H_lost
			x[1][:, :, :, cripple] = 0 # H_gain
			# Set atom attributes within bond attributes to zero
			x[2][:, :, :, cripple] = 0 # Bond_lost
			x[2][:, :, :, cripple + offset] = 0 # Bond_lost
			x[3][:, :, :, cripple] = 0 # Bond_gain
			x[3][:, :, :, cripple + offset] = 0 # Bond_gain
			print('Set index {} and {} to zero'.format(cripple, cripple + offset))

		preds = model.predict(x, batch_size = 20)

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
	parser.add_argument('--tag', type = str, default = int(time.time()),
						help = 'Tag for this model')
	parser.add_argument('--fold', type = int, default = 5, 
						help = 'Which fold of the 5-fold CV is this? Defaults 5')
	parser.add_argument('data', type = str,
		                help = 'Data folder with data files')
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
	N_c = 500 # number of candidate edit sets
	N_e = 5 # maximum number of edits per class

	THIS_FOLD_OUT_OF_FIVE = int(args.fold)
	tag = args.tag + ' fold{}'.format(args.fold)

	print('Reloading from file')
	model = model_from_json(open(os.path.join(FROOT, 'model{}.json'.format(tag))).read())
	model.compile(loss = 'categorical_crossentropy', 
		optimizer = SGD(lr = 0.01, decay = 1e-4, momentum = 0.9),
		metrics = ['accuracy']
	)
	model.load_weights(os.path.join(FROOT, 'weights{}.h5'.format(tag)))
	
	fid = open(os.path.join(FROOT, 'input_masking_{}.csv'.format(tag)), 'w')
	fid.write('crippled_index\ttrain_acc\tval_acc\n')

	for i in range(F_atom):
                print('CRIPPLING {}'.format(i))
                # Cripple index i
		(train_acc, val_acc) = get_accuracy(model, x_files, y_files, z_files, tag = tag, split_ratio = 0.8, cripple = i)
		fid.write('{}\t{}\t{}\n'.format(i, train_acc, val_acc))
