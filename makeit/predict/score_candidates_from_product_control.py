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
import rdkit.Chem.AllChem as AllChem
import theano.tensor as T
from scipy.sparse import coo_matrix
import cPickle as pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt    # for visualization
import scipy.stats as ss
import itertools



def build(F_atom = 1, F_bond = 1, N_e = 5, N_h1 = 100, N_h2 = 50, N_h3 = 0, N_c = 500, inner_act = 'tanh',
		l2v = 0.01, lr = 0.0003, N_hf = 0):
	'''
	Builds the feed forward model.

	N_c:  maximum number of candidates for each sample (padded length)
	inner_act: activation function 
	'''

	net_sum = Input(shape = (N_c, 1024))

	if N_hf > 0:
		feature_to_feature = Dense(N_hf, activation = inner_act, W_regularizer = l2(l2v))
		net_sum_h = TimeDistributed(feature_to_feature)(net_sum)
	else:
		net_sum_h = net_sum 

	feature_to_score = Dense(1, activation = 'linear', W_regularizer = l2(l2v))
	unscaled_score = TimeDistributed(feature_to_score)(net_sum_h)

	unscaled_score_flat = Flatten()(unscaled_score)

	score = Activation('softmax')(unscaled_score_flat)
	
	model = Model(input = [net_sum], output = [score])

	# Now compile
	sgd = SGD(lr = lr, decay = 1e-4, momentum = 0.9)
	# model.compile(loss = 'binary_crossentropy', optimizer = sgd, 
	# 	metrics = ['binary_accuracy'])
	model.compile(loss = 'categorical_crossentropy', optimizer = sgd, 
		metrics = ['accuracy'])

	model.summary()

	return model

def train(model, x_files, y_files, z_files, tag = '', split_ratio = 0.8):

	hist_fid = open(os.path.join(FROOT, 'hist{}.csv'.format(tag)), 'a')
	hist_fid.write('epoch,filenum,loss,val_loss,acc,val_acc\n')
	try:
		for epoch in range(nb_epoch):
			print('>>> EPOCH {}/{} <<<'.format(epoch + 1, nb_epoch))
			average_loss = 0; average_acc = 0;
			average_val_loss = 0; average_val_acc = 0;
			# Get each set of data
			for fnum in range(len(x_files)):
				print('    running file {}/{}'.format(fnum+1, len(x_files)))
				
				if len(x_files) > 1 or epoch == 0: # get new data
					z = get_z_data(z_files[fnum])
					x = get_x_data(x_files[fnum], z)
					y = get_y_data(y_files[fnum])
					z = rearrange_for_5fold_cv(z)
					
				hist = model.fit(x, y, 
					nb_epoch = 1, 
					batch_size = batch_size, 
					validation_split = (1 - split_ratio),
					verbose = 0)

				hist_fid.write('{},{},{},{},{},{}\n'.format(
					epoch + 1, fnum, hist.history['loss'][0], hist.history['val_loss'][0],
					hist.history['acc'][0], hist.history['val_acc'][0]
				))
				# print('loss: {}, val_loss: {}, acc: {}, val_acc: {}'.format(
				#	hist.history['loss'][0], hist.history['val_loss'][0],
				# 	hist.history['acc'][0], hist.history['val_acc'][0]
				#))
				average_loss += hist.history['loss'][0]
				average_acc += hist.history['acc'][0]
				average_val_loss += hist.history['val_loss'][0]
				average_val_acc += hist.history['val_acc'][0]
			print('    loss:     {:8.4f}, acc:     {:5.4f}'.format(average_loss/len(x_files), average_acc/len(x_files)))
			print('    val_loss: {:8.4f}, val_acc: {:5.4f}'.format(average_val_loss/len(x_files), average_val_acc/len(x_files)))
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
	print('REARRANGING FOR CV FOLD {}'.format(THIS_FOLD_OUT_OF_FIVE))
	chunks = [lst[i:i+N] for i in range(0, len(lst), N)]
	new_lst = list(itertools.chain.from_iterable(chunks[:(THIS_FOLD_OUT_OF_FIVE-1)])) + \
		list(itertools.chain.from_iterable(chunks[THIS_FOLD_OUT_OF_FIVE:])) + \
		list(itertools.chain.from_iterable(chunks[(THIS_FOLD_OUT_OF_FIVE-1):THIS_FOLD_OUT_OF_FIVE]))
	if type(lst) == type(np.array(1)):
		return np.array(new_lst)
	return new_lst


def get_x_data(fpath, z):
	'''
	Reads the candidate SMILES and returns a list of fingerprints

	for the control experiment
	'''

	# Check if we've already populated these matrices once
	if os.path.isfile(fpath + '_FP_processed'):
		with open(fpath + '_FP_processed', 'rb') as infile:
			FPs = pickle.load(infile)
			return rearrange_for_5fold_cv(FPs)

	# Otherwise, we need to do it...
	with open(fpath, 'rb') as infile:
		all_candidate_smiles = pickle.load(infile)
		N = len(z)
		FPs = np.zeros((N, N_c, 1024))

		for i, candidates in enumerate(all_candidate_smiles):
			for j, candidate in enumerate(candidates):
				prod = Chem.MolFromSmiles(str(candidate))
				fp = np.array(AllChem.GetMorganFingerprintAsBitVect(prod, 2, nBits = 1024))
				FPs[i, j, :] = fp

		# Dump file so we don't have to do that again
		with open(fpath + '_FP_processed', 'wb') as outfile:
			pickle.dump(FPs, outfile, pickle.HIGHEST_PROTOCOL)
		print('Converted {} to FPs for the first (and only) time'.format(fpath))
		return rearrange_for_5fold_cv(FPs)

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

def pred_histogram(model, x_files, y_files, z_files, tag = '', split_ratio = 0.8):
	'''
	Given a trained model and a list of samples, this function creates a 
	histogram of the pseudo-probabilities assigned to the true reaction 
	examples (i.e., the correct product)
	'''

	fid = open(os.path.join(FROOT, 'probs{}.dat'.format(tag)), 'w')
	fid.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
		'reaction_smiles', 'train/val', 
		'true_edit', 'prob_true_edit', 
		'predicted_edit', 'prob_predicted_edit',
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
		x = get_x_data(x_files[fnum], z)
		y = get_y_data(y_files[fnum])
		z = rearrange_for_5fold_cv(z)

		# Getting unprocessed x data
		with open(x_files[fnum], 'rb') as infile:
			all_edits = pickle.load(infile)

		preds = model.predict(x, batch_size = batch_size)
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
			most_likely_edit_i = np.argmax(pred)
			if most_likely_edit_i >= len(edits): # no reaction?
				most_likely_edit = 'no_reaction'
			else:
				most_likely_edit = edits[most_likely_edit_i]
			fid.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
				z[i], dataset, 
				edits[np.argmax(y[i,:])], trueprob, 
				most_likely_edit, np.max(pred),
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
	args = parser.parse_args()

	x_files = sorted([os.path.join(args.data, dfile) \
					for dfile in os.listdir(args.data) \
					if 'candidate_smiles' in dfile and '_processed' not in dfile])
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
	N_c = int(args.Nc) # number of candidate edit sets
	N_e = 5 # maximum number of edits per class

	THIS_FOLD_OUT_OF_FIVE = int(args.fold)
	tag = args.tag + ' fold{}'.format(args.fold)

	if bool(args.retrain):
		print('Reloading from file')
		model = model_from_json(open(os.path.join(FROOT, 'model{}.json'.format(tag))).read())
		model.compile(loss = 'categorical_crossentropy', 
			optimizer = SGD(lr = lr, decay = 1e-4, momentum = 0.9),
			metrics = ['accuracy']
		)
		model.load_weights(os.path.join(FROOT, 'weights{}.h5'.format(tag)))
	else:
		model = build(F_atom = F_atom, F_bond = F_bond, N_e = N_e, N_c = N_c, N_h1 = N_h1, N_h2 = N_h2, N_h3 = N_h3, N_hf = N_hf, l2v = l2v, lr = lr)
		with open(os.path.join(FROOT, 'model{}.json'.format(tag)), 'w') as outfile:
        	        outfile.write(model.to_json())

	if bool(args.test):
		pred_histogram(model, x_files, y_files, z_files, tag = tag, split_ratio = 0.8)
		quit(1)
	elif bool(args.visualize):
		visualize_weights(model, tag)
		quit(1)

	hist = train(model, x_files, y_files, z_files, tag = tag, split_ratio = 0.8)
	model.save_weights(os.path.join(FROOT, 'weights{}.h5'.format(tag)), overwrite = True) 

	pred_histogram(model, x_files, y_files, z_files, tag = tag, split_ratio = 0.8)
	#visualize_weights(model, tag)
