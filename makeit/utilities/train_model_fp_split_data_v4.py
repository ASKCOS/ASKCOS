
#import files
from __future__ import division

import theano
from theano.tensor import lt,le,eq,gt,ge

import numpy as np#something
import datetime
import time
import os
import sys
import argparse

import h5py # needed for save_weights, fails otherwise

from keras import backend as K 
from keras.models import Sequential, Model, model_from_json
from keras.layers import Dense, Activation, Input
from keras.layers.core import Flatten, Permute, Reshape, Dropout, Lambda, RepeatVector
from keras.layers.wrappers import TimeDistributed
from keras.layers import Merge, merge, activations
from keras.layers.merge import Dot, Add
from keras.optimizers import SGD, Adam, Adadelta
from keras.layers.convolutional import Convolution1D, Convolution2D
from keras.regularizers import l2
from keras.utils.np_utils import to_categorical
from keras.utils.generic_utils import func_dump
from keras.utils.generic_utils import func_load
#import keras.engine.topology 
from keras.engine.topology import Layer
import pickle
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import os
from sklearn.metrics import roc_auc_score
from rdkit import RDLogger
from scipy import sparse
import pandas as pd
from tqdm import tqdm
import random

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
####utilities
def set_keras_backend(backend):

    if K.backend() != backend:
        os.environ['KERAS_BACKEND'] = backend
        reload(K)
        assert K.backend() == backend

# def load_sparse_csr(filename):
#     loader = np.load(filename)
#     return sparse.csr_matrix((  loader['data'], loader['indices'], loader['indptr']),
#                          shape = loader['shape'])
#define a highway layer
#######################tried to embed DENSE layer, failed...############################
# class Highway(Layer):

# 	def __init__(self, activation, l2v = 0.01, **kwargs):
# 		super(Highway, self).__init__(**kwargs)
# 		self.activation = activation
# 		self.transform_actv = activations.get('sigmoid')
# 		self.l2v = l2v

# 	def build(self, input_shape):
# 		#weights of the dense layer

# 		self.kernel_T = self.add_weight(name = 'kernel_T',
# 								 shape = (input_shape[1],input_shape[1]),
# 								 initializer ='glorot_uniform',
# 								 trainable = True)
# 		self.bias_T = self.add_weight(name = 'bias_T',
# 								 shape = (input_shape[1],),
# 								 initializer ='zeros',
# 								 trainable = True)
# 		self.input_dim = input_shape[1]
# 		# print(self.input_dim)
# 		super(Highway, self).build(input_shape)
	
# 	def call(self, x):
# 		transform_fun = Dense(self.input_dim, activation = self.activation,kernel_regularizer = l2(self.l2v))(x)
# 		transform_gate = self.transform_actv(K.bias_add(K.dot(x,self.kernel_T), self.bias_T))
# 		carry_gate = K.ones(self.input_dim,) - transform_gate
# 		output = transform_fun*transform_gate + x*carry_gate
# 		return output
#############################################################################

class Highway_self(Layer):

	def __init__(self, activation = 'elu',**kwargs):
		super(Highway_self, self).__init__(**kwargs)
		self.activation = activations.get(activation)
		self.transform_actv = activations.get('sigmoid')
		

	def build(self, input_shape):
		#weights of the dense layer
		self.kernel = self.add_weight(name = 'kernel',
								 shape = (input_shape[1],input_shape[1]),
								 initializer ='glorot_uniform',
								 trainable = True)
		self.bias = self.add_weight(name = 'bias',
								 shape = (input_shape[1],),
								 initializer ='zeros',
								 trainable = True)
		self.kernel_T = self.add_weight(name = 'kernel_T',
								 shape = (input_shape[1],input_shape[1]),
								 initializer ='glorot_uniform',
								 trainable = True)
		self.bias_T = self.add_weight(name = 'bias_T',
								 shape = (input_shape[1],),
								 initializer ='zeros',
								 trainable = True)
		self.input_dim = input_shape[1]
		# print(self.input_dim)
		super(Highway_self, self).build(input_shape)
	
	def call(self, x):
		transform_fun = self.activation(K.bias_add(K.dot(x,self.kernel), self.bias))
		transform_gate = self.transform_actv(K.bias_add(K.dot(x,self.kernel_T), self.bias_T))
		carry_gate = K.ones(self.input_dim,) - transform_gate
		output = transform_fun*transform_gate + x*carry_gate
		return output

	def compute_output_shape(self, input_shape):
		return (input_shape[0],input_shape[1])

# 	def compute_output_shape(self, input_shape):
# 		return (input_shape[0],input_shape[1])

#function for creating finger prints fingerprints
def create_rxn_Morgan2FP(rsmi, psmi, rxnfpsize = 2048, pfpsize=2048, useFeatures = False,calculate_rfp = True):
    """Create a rxn Morgan (r=2) fingerprint as bit vector from SMILES string lists of reactants and products"""
    # Modified from Schneider's code (2014)
    if calculate_rfp is True:
	    rsmi = rsmi.encode('utf-8')
	    try:
	    	mol = Chem.MolFromSmiles(rsmi)
	    except Exception as e:
	    	
	    	return
	    # if mol is None:
	    # 	print(react)
	    try:
	        fp_bit = AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits = rxnfpsize, useFeatures=False)
	        fp = np.empty(rxnfpsize,dtype = 'int8')
	        DataStructs.ConvertToNumpyArray(fp_bit,fp)
	        # print(fp.dtype)
	        # fp = np.asarray(fp_bit)
	        # fp = AllChem.GetMorganFingerprint(mol=mol, radius=2, useFeatures=useFeatures)

	    except Exception as e:
	        print("Cannot build reactant fp due to {}".format(e))

	        return
	        
	    rfp = fp
    else:
	    rfp = None

    psmi = psmi.encode('utf-8')
    try:
    	mol = Chem.MolFromSmiles(psmi)
    except Exception as e:
    	print(psmi)
    	return
    # if mol is None:
    # 	print(product)
    try:
        fp_bit = AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits = pfpsize, useFeatures=False)
        fp = np.empty(pfpsize,dtype = 'int8')
        DataStructs.ConvertToNumpyArray(fp_bit,fp)
        # fp = np.asarray(fp_bit)
        # fp = AllChem.GetMorganFingerprint(mol=mol, radius=2, useFeatures=useFeatures)

    except Exception as e:
    	print("Cannot build product fp due to {}".format(e))
    	return
        
    pfp = fp
    # pfp_for_rxn = pfp
    # for product in psmi:
    # 	product = product.encode('utf-8')
    #     mol = Chem.MolFromSmiles(product)
    #     if mol is None:
    #     	print(product)
    #     try:
    #         fp = np.array(
    #             AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits=rxnfpsize, useFeatures=useFeatures))
    #     except Exception as e:
    #         print("Cannot build product fp due to {}".format(e))
    #     if pfp_for_rxn is None:
    #         pfp_for_rxn = fp
    #     else:
    #         pfp_for_rxn += fp
    # if pfp_for_rxn is not None and rfp is not None:
    #     rxnfp = pfp_for_rxn - rfp
    return [pfp, rfp]


#function for loading and partitioning data
def load_and_partition_data(pfp_csr_matrix, rfp_csr_matrix, outcomes_list, pfp_label_list, split_ratio, batch_size):
	N_samples = len(outcomes_list)
	N_train = int(N_samples * split_ratio[0])
	N_val	= int(N_samples * split_ratio[1])
	N_test  = N_samples - N_train - N_val
	print('Total number of samples: {}'.format(N_samples))
	print('Training   on {}% - {}'.format(split_ratio[0]*100, N_train))
	print('Validating on {}% - {}'.format(split_ratio[1]*100, N_val))
	print('Testing    on {}% - {}'.format((1-split_ratio[1]-split_ratio[0])*100, N_test))


	return {
		'N_samples': N_samples,
		'N_train': N_train,
		#
		'train_generator': batch_data_generator(pfp_csr_matrix, rfp_csr_matrix, outcomes_list, 0, N_train, batch_size),
		'train_label_generator': batch_label_generator(pfp_label_list, outcomes_list, 0, N_train, batch_size),
		'train_nb_samples': N_train,
		#
		'val_generator': batch_data_generator(pfp_csr_matrix, rfp_csr_matrix, outcomes_list, N_train, N_train + N_val, batch_size),
		'val_label_generator': batch_label_generator(pfp_label_list, outcomes_list, N_train, N_train + N_val, batch_size),
		'val_nb_samples': N_val,
		#
		'test_generator': batch_data_generator(pfp_csr_matrix, rfp_csr_matrix, outcomes_list, N_train + N_val, N_samples, batch_size),
		'test_label_generator': batch_label_generator(pfp_label_list, outcomes_list, N_train + N_val, N_samples, batch_size),
		'test_nb_samples': N_test,
		#
		#
		'batch_size': batch_size,
	}


#batch data generator
def batch_data_generator( pfp_csr_matrix, rfp_csr_matrix, outcomes_list, start_at, end_at, batch_size):
	
	#X_train_batch = []
	# y_train_batch = []


	while True:
		for start_index in range(start_at, end_at,batch_size):
			end_index = min(start_index + batch_size, end_at)
		##anotherway of starting the loop
		# start_index = start_at-1
		# while start_index < end_at-1:

			# X_train_batch[start_index:end_index] = \
			y_train_batch = outcomes_list[start_index:end_index]
			pfp_matrix = pfp_csr_matrix[start_index:end_index,:].todense()
			rfp_matrix = rfp_csr_matrix[start_index:end_index,:].todense()
			
			# print(pfp_train_batch.shape)
			# print(rxnfp_train_batch.shape)
			# print(y_train_batch.shape)
			# print(type(X_train_batch))
			# # print(X_train_batch.shape)
			# y_train_batch = \
			# 	np.asarray([x[2] for x in data_file[start_index:end_index]], dtype='float32')
			pfp_train_batch = np.asarray(pfp_matrix,dtype = 'float32')
			rfp_train_batch = np.asarray(rfp_matrix,dtype = 'float32')

			rxnfp_train_batch = pfp_train_batch - rfp_train_batch
			y_train_batch = np.asarray(y_train_batch,dtype = 'float32')
			# print(y_train_batch.shape)
			# print(type(y_train_batch))
			# print(y_train_batch.shape)
			## Reshape happening here?
			# X_train_batch[start_index:end_index] = [a.astype('float32') for a in X_train_batch[start_index:end_index]]
			# y_train_batch[start_index:end_index] = [a.astype('float32') for a in y_train_batch[start_index:end_index]]
			# print(pfp_train_batch.shape)
			# print(rxnfp_train_batch.shape)
			yield ([pfp_train_batch, rxnfp_train_batch],y_train_batch)


def batch_label_generator(pfp_label_list, outcomes_list, start_at, end_at, batch_size):
	# while True:
		# for start_index in range(start_at, end_at,batch_size):
		# 	end_index = min(start_index + batch_size, end_at)

			# rxn_smiles = [] 
			# rxn_true = []

			# for x in data_file[start_index:end_index]:
			# 	rxn_smiles.append('>>'.join([x[0],x[1]]))
			# 	rxn_true.append(x[2])
	while True:
		for start_index in range(start_at, end_at,batch_size):
			end_index = min(start_index + batch_size, end_at)

		##anotherway of starting the loop
		# start_index = start_at-1
		# while start_index < end_at-1:

			# X_train_batch[start_index:end_index] = \
			rxn_id = []
			rxn_true = []
			for i in range(start_index,end_index):
				rxn_id.append(pfp_label_list[i])
				rxn_true.append(outcomes_list[i])
				
	# while True:
	# 	# for start_index in range(start_at, end_at,batch_size):
	# 	# 	end_index = min(start_index + batch_size, end_at)
	# 	##anotherway of starting the loop
	# 	start_index = start_at-1
	# 	while start_index < end_at-1:
	# 		rxn_smiles = []
	# 		rxn_true = []
	# 		while len(rxn_smiles)<batch_size:
	# 			start_index += 1 
	# 			x = data_file[start_index]
	# 			try:
	# 				[pfp, rxnfp] = create_rxn_Morgan2FP(x[0].split('.'),x[1].split('.'))
	# 				rxn_smiles.append('>>'.join([x[0],x[1]]))
	# 				rxn_true.append(x[2])
	# 			except Exception as e:
	# 				continue

			yield (rxn_id, rxn_true)

def multiple_batch_data_generator(data_generator_list):
	while True:
		pfp_train_batch = []
		rxnfp_train_batch = []
		y_train_batch = []
		for i in range(len(data_generator_list)):
			(x,y) = data_generator_list[i].next()
			if pfp_train_batch == []:
				pfp_train_batch = x[0]
			else:
				pfp_train_batch = np.append(pfp_train_batch,x[0],0)
			if rxnfp_train_batch == []:
				rxnfp_train_batch = x[1]
			else:
				rxnfp_train_batch = np.append(rxnfp_train_batch,x[1],0)
			if y_train_batch == []:
				y_train_batch = y
			else:
				y_train_batch = np.append(y_train_batch,y)
			# pfp_train_batch = np.asarray(pfp_train_batch,dtype = 'float32')
			# rxnfp_train_batch = np.asarray(rxnfp_train_batch,dtype = 'float32')
			# y_train_batch = np.asarray(y_train_batch,dtype = 'float32')
		yield ([pfp_train_batch,rxnfp_train_batch],y_train_batch)

def multiple_batch_label_generator(label_generator_list):
	while True:
		batch_label = [[],[]]
		for i in range(len(label_generator_list)):
			(x,y) = label_generator_list[i].next()
			batch_label[0] = np.append(batch_label[0],x,0)
			batch_label[1] = np.append(batch_label[1],y,0)
		yield batch_label
#build model structure

def pos_ct(y_true, y_pred):
	pos_pred = K.sum(gt((K.clip(y_pred, 0, 1)),0.5))
	return pos_pred
def true_pos(y_true, y_pred):
	true_pos_ct = K.sum(gt((K.clip(y_pred*y_true, 0, 1)),0.5))
	return true_pos_ct
def real_pos(y_true, y_pred):
	real_pos_ct = K.sum(gt((K.clip(y_true, 0, 1)),0.5))
	return real_pos_ct

def build(pfp_len = 2048, rxnfp_len = 2048,l2v = 0.01):
	input_pfp = Input(shape = (pfp_len,))
	input_rxnfp = Input(shape = (rxnfp_len,))
	
	input_pfp_h1 = Dense(1024, activation = 'elu')(input_pfp)
	input_pfp_h2 = Dropout(0.3)(input_pfp_h1)
	input_pfp_h3 = Highway_self(activation = 'elu')(input_pfp_h2)
	input_pfp_h4 = Highway_self(activation = 'elu')(input_pfp_h3)
	input_pfp_h5 = Highway_self(activation = 'elu')(input_pfp_h4)
	input_pfp_h6 = Highway_self(activation = 'elu')(input_pfp_h5)
	input_pfp_h7 = Highway_self(activation = 'elu')(input_pfp_h6)

	# input_pfp_h3 = Dense(1024,activation = 'elu')(input_pfp_h2)
	# input_pfp_h4 = Dense(1024,activation = 'elu')(input_pfp_h3)
	# input_pfp_h5 = Dense(1024,activation = 'elu')(input_pfp_h4)
	# input_pfp_h6 = Dense(1024,activation = 'elu')(input_pfp_h5)
	# input_pfp_h7 = Dense(1024,activation = 'elu')(input_pfp_h6)
	input_rxnfp_h1 = Dense(1024, activation = 'elu')(input_rxnfp)
	merged_h1 = Dot(axes = 1, normalize=False)([input_pfp_h7,input_rxnfp_h1])
	
	# scale = Dense(1, activation = 'tanh')(merged_h1)
#	model.add(Flatten())
	
	output= Dense(1, activation = 'sigmoid')(merged_h1)
	model = Model([input_pfp,input_rxnfp],output)
	#print(model.output_shape)

	model.count_params()
	model.summary()
	
	adam = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.0)

	# 	# return K.sum(y_true)
	# def false_pos_ct(y_true, y_pred):
	# 	return K.sum(gt(y_pred,0.5) and eq(y_true,0))
	# def false_neg_ct(y_true, y_pred):
	# 	return K.sum(lt(y_pred,0.5) and eq(y_true,1))

	model.compile(loss='binary_crossentropy', optimizer='adam', metrics = ['acc', pos_ct, true_pos, real_pos])
	return model

#train model
def train_bin_class(model, pfp_csr_matrices, rfp_csr_matrices, outcomes_lists, pfp_label_lists, split_ratio, class_weight, batch_size):
	train_generator_list = []
	val_generator_list = []
	train_nb_samples = 0
	val_nb_samples = 0
	data = [0]*len(outcomes_lists)
	print(len(data))
	batch_size_total = 0
	for i in range(len(outcomes_lists)):
		data[i] = load_and_partition_data(pfp_csr_matrices[i], rfp_csr_matrices[i], outcomes_lists[i], pfp_label_lists[i], split_ratio, batch_size[i])
		train_generator_list.append(data[i]['train_generator'])
		val_generator_list.append(data[i]['val_generator'])
		train_nb_samples += data[i]['train_nb_samples']
		val_nb_samples += data[i]['val_nb_samples']
		batch_size_total += data[i]['batch_size']
	train_generator = multiple_batch_data_generator(train_generator_list)
	val_generator = multiple_batch_data_generator(val_generator_list)
	
	(x,y) = train_generator.next()
	# print(x[0].shape,x[1].shape,y.shape)
	# print(train_generator.next())
	# print(val_generator.next())
	# wait = raw_input("press to cont")
	# Add additional callbacks
	from keras.callbacks import ModelCheckpoint, CSVLogger, EarlyStopping, ReduceLROnPlateau
	reduce_lr = ReduceLROnPlateau(monitor='loss', factor=0.2,
              patience=2, min_lr=0.0001)
	callbacks = [
		# ModelCheckpoint(WEIGHTS_FPATH), # save every epoch
		EarlyStopping(patience = 5),
		CSVLogger('training.log'),
		reduce_lr
		]

	try:
		hist = model.fit_generator(train_generator, 
			verbose = 1,
			validation_data = val_generator,
			steps_per_epoch = np.ceil(train_nb_samples/batch_size_total),
			epochs = nb_epoch, 
			callbacks = callbacks,
			validation_steps = np.ceil(val_nb_samples/batch_size_total),
			class_weight = class_weight,
		)

	except KeyboardInterrupt:
		print('Stopped training early!')

#test model
# def test_bin_class(model, pfp_csr_matrix,pfp_label_list, outcomes_list, rfp_df, split_ratio):
# 	data = load_and_partition_data(pfp_csr_matrix,pfp_label_list, outcomes_list, rfp_df, split_ratio)
# 	TEST_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/results/test_results.dat"
# 	SENS_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/results/sens_to_thres.dat"
# 	fid = open(TEST_FPATH, 'w')
# 	sens_to_thres =open(SENS_FPATH, 'w')
# 	def test_on_set(fid, dataset, data_generator, label_generator, num_batches):
# 		print('Testing on {} data'.format(dataset))

# 		outcome_true = []
# 		outcome_pred = []
# 		corr = 0

# 		#test for different thresholds
# 		true_pos = [0,0,0,0,0]
# 		true_neg = [0,0,0,0,0]
# 		false_pos = [0,0,0,0,0]
# 		false_neg = [0,0,0,0,0]

# 		for batch_num in range(num_batches):
# 			(x, y_true) = data_generator.next()
# 			(rxn_id, rxn_true) = label_generator.next()

# 			y_pred = model.predict_on_batch(x)

# 			for i in range(y_pred.shape[0]):
# 				outcome_true.append(y_true[i])
# 				outcome_pred.append(y_pred[i])
# 				# print(y_true[i],np.rint(y_pred[i]))
# 				# print corr
# 				if (y_true[i]- np.rint(y_pred[i]) == 0.0):
# 					corr+=1

# 				if (y_true[i] + np.rint(y_pred[i])>=1):
# 					fid.write('{}\t{}\t{}\t{}\t{}\n'.format(
# 								dataset,
# 								rxn_id[i],
# 								y_true[i],
# 								rxn_true[i],
# 								y_pred[i],
# 							))
# 		class_thresholds = [0.1,0.2,0.3,0.4,0.5]
# 		sens_to_thres.write('dataset\tthreshold\ttrue_pos\tfalse_neg\tfalse_pos\ttrue_neg\n')

# 		for i in range(5):
# 			for j in range(len(outcome_true)):
# 				if (outcome_true[j] == 1 and outcome_pred[j]>class_thresholds[i]):
# 					true_pos[i] +=1
# 				elif (outcome_true[j] == 1 and outcome_pred[j]<=class_thresholds[i]):
# 					false_neg[i]+=1
# 				elif (outcome_true[j] == 0 and outcome_pred[j]<=class_thresholds[i]):
# 					true_neg[i]+=1
# 				elif (outcome_true[j] == 0 and outcome_pred[j]>class_thresholds[i]):
# 					false_pos[i]+=1
# 			sens_to_thres.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
# 								dataset,
# 								class_thresholds[i],
# 								true_pos[i],
# 								false_neg[i],
# 								false_pos[i],
# 								true_neg[i]
# 								))


# 		# print(corr, len(outcome_true))
# 		accuracy = corr/len(outcome_true)
# 		# print(accuracy)
# 		# print(outcome_true[:10])
# 		# print(outcome_pred[:10])
# 		roc_score = roc_auc_score(outcome_true, outcome_pred)


# 		fid.write('overall: {}\t{}\t{}\n\n\n\n'.format(
# 					dataset, 
# 					accuracy,
# 					roc_score,
# 				))


# 		return accuracy, roc_score

# 	train_accu, train_auc = test_on_set(fid, 'train', data['train_generator'], data['train_label_generator'], 
# 		int(np.ceil(data['train_nb_samples']/float(data['batch_size'])))
# 	)
# 	val_accu, val_auc = test_on_set(fid, 'val', data['val_generator'], data['val_label_generator'], 
# 		int(np.ceil(data['val_nb_samples']/float(data['batch_size'])))
# 	)
# 	test_accu, test_auc = test_on_set(fid, 'test', data['test_generator'], data['test_label_generator'], 
# 		int(np.ceil(data['test_nb_samples']/float(data['batch_size'])))
# 	)

# 	fid.close()
# 	sens_to_thres.close()

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('--nb_epoch', type = int, default = 100,
						help = 'Number of epochs to train for, default 100')
	parser.add_argument('--batch_size', type = int, default = 20,
						help = 'Batch size, default 20')
	parser.add_argument('--retrain', type = bool, default = False,
						help = 'whether to train the model, default False')
	args = parser.parse_args()

	nb_epoch           = int(args.nb_epoch)
	batch_size_input   = int(args.batch_size)
	retrain			   = bool(args.retrain)

	WEIGHTS_FPATH = "model_weights1.h5"
	SPLIT_RATIO = (0.8, 0.1)
	# nb_epoch = 1
	# batch_size = 32
	set_keras_backend("theano")
	model = build(pfp_len = 2048, rxnfp_len = 2048, l2v = 1e-9)
	
	print("model built")
	# data_file_raw = "bin_class_dataset_medium.pickle"

	# with open(data_file_raw,"r") as f:
	#  	data_file = pickle.load(f)
	# print(len(data_file))
	rid_lists = []
	outcomes_lists = []
	pfp_csr_matrices = []
	rfp_csr_matrices = []
	PFP1_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/reaxys/pfp_dataset_1.npz"
	RFP1_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/reaxys/rfp_dataset_1.npz"
	RID1_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/reaxys/rid_1.pickle"
	BIN1_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/reaxys/outcome_dataset_1.pickle"

	with open(RID1_FPATH,"r") as RLB:
		rid_lists.append(pickle.load(RLB))
	outcomes_lists.append([1.0]*len(rid_lists[0]))
	# with open(BIN1_FPATH,"r") as BIN:
	# 	outcomes_lists.append(pickle.load(BIN)) 
	pfp_csr_matrices.append(sparse.load_npz(PFP1_FPATH))
	rfp_csr_matrices.append(sparse.load_npz(RFP1_FPATH))

	PFP2_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_1_200k_prediction/pfp_dataset_1.npz"
	RFP2_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_1_200k_prediction/rfp_dataset_1.npz"
	RID2_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_1_200k_prediction/rid_1.pickle"
	BIN2_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_1_200k_prediction/outcome_dataset_1.pickle"

	with open(RID2_FPATH,"r") as RLB:
		rid_lists.append(pickle.load(RLB))
	with open(BIN2_FPATH,"r") as BIN:
		outcomes_lists.append(pickle.load(BIN)) 
	pfp_csr_matrices.append(sparse.load_npz(PFP2_FPATH))
	rfp_csr_matrices.append(sparse.load_npz(RFP2_FPATH))

	PFP3_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_2_100k_prediction/pfp_dataset_2.npz"
	RFP3_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_2_100k_prediction/rfp_dataset_2.npz"
	RID3_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_2_100k_prediction/rid_2.pickle"
	BIN3_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_2_100k_prediction/outcome_dataset_2.pickle"

	with open(RID3_FPATH,"r") as RLB:
		rid_lists.append(pickle.load(RLB))
	with open(BIN3_FPATH,"r") as BIN:
		outcomes_lists.append(pickle.load(BIN)) 
	pfp_csr_matrices.append(sparse.load_npz(PFP3_FPATH))
	rfp_csr_matrices.append(sparse.load_npz(RFP3_FPATH))

	print("raw data loaded, converting rfp matrix to dense...\n")

	nb_sample_1 = len(outcomes_lists[0])
	nb_sample_2 = len(outcomes_lists[1])
	nb_sample_3 = len(outcomes_lists[2])
	batch_size = [0]*3
	batch_size[0] = int(float(nb_sample_1)/(nb_sample_1 + nb_sample_2 + nb_sample_3)*batch_size_input)
	batch_size[1] = int(float(nb_sample_2)/(nb_sample_1 + nb_sample_2 + nb_sample_3)*batch_size_input)
	batch_size[2] = batch_size_input - batch_size[0] - batch_size[1]
	# batch_size = [batch_size_input]
	print(batch_size)
	print(len(outcomes_lists),len(pfp_csr_matrices),len(rfp_csr_matrices),len(rid_lists))
	# pfp_label_list = pfp_label_list[0:70000000]
	# print("label_ok")
	# outcomes_list = outcomes_list[0:70000000]
	# print("outcome ok..")
	# pfp_csr_matrix = pfp_csr_matrix[0:700000,:]
	# print("csr_matrix_ok...")



	print("raw data loaded, converting rfp matrix to dense...\n")
	print(datetime.datetime.now())

	# pfp_matrix = pfp_csr_matrix.todense()

	#only convert rfp_matrix to dense



	# print('merge complete...')
	# df_shuffled['rxnfp'] = df_shuffled['pfp']-df_shuffled['rfp']
	# # print(df_shuffled['rxnfp'].head(5))
	# print(df_shuffled['pfp'].head(5))
	# pfp_label_list = list(df_shuffled['rxn_id'])
	# pfp_list = list(df_shuffled['pfp'])
	# outcomes_list = list(df_shuffled['outcome'])
	# rfp_list = list(df_shuffled['rxnfp'])
	# print('shuffling complete...')
	# print(pfp_list[0:20])
	# print(pfp_label_list[0:20])
	# print(outcomes_list[0:20])
	# # for test purpose use small subset of data
	# # data_file = data_file[:1000]
	# # print(len(data_file))

	
	# train model if retrain is true
	wt_file = h5py.File(WEIGHTS_FPATH,'w')

	# PREV_WEIGHTS_PATH = "previous_model_weights/model_weights1.h5"
	# prev_wt_file = h5py.File(PREV_WEIGHTS_PATH,'r')

	class_weight = None
	# class_weight = {0 : 1.,
	# 1: 100.}
	retrain = True
	# if retrain is True:
	# weights = []
	# for i in range(len(prev_wt_file.keys())):
	#     weights.append(prev_wt_file['weight'+str(i)][:])
	# model.set_weights(weights)

	train_bin_class(model, pfp_csr_matrices, rfp_csr_matrices, outcomes_lists, rid_lists, SPLIT_RATIO, class_weight, batch_size)
	model.save('my_model.h5')

	from keras.models import load_model
	model = load_model('my_model.h5', custom_objects = {'Highway_self':Highway_self,'pos_ct':pos_ct, 'true_pos':true_pos, 'real_pos':real_pos})

	# weights = model.get_weights()
	# for i in range(len(weights)):
	# 	wt_file.create_dataset('weight'+str(i),data=weights[i])
	# # else:
	# # 	print('Directly loading weights and make predictions...')

	# weights = []
	# for i in range(len(wt_file.keys())):
	#     weights.append(wt_file['weight'+str(i)][:])
	# model.set_weights(weights)

	# wt_file.close()
	# print('weights loaded...')
	#model.save_weights(WEIGHTS_FPATH, overwrite = True) #doesnt work for bidirectional network
	# test_bin_class(model, pfp_csr_matrix,pfp_label_list, outcomes_list, rfp_df, SPLIT_RATIO)


	################
	# Code for retriving weights
	# file=h5py.File(WEIGHTS_FPATH,'r')
	# weight = []
	# for i in range(len(file.keys())):
	#     weight.append(file['weight'+str(i)][:])
	# model.set_weights(weight)
	# ... ...
	# ... ...
	# file = h5py.File(fileName,'w')
	# weight = model.get_weights()
	# for i in range(len(weight)):
	#     file.create_dataset('weight'+str(i),data=weight[i])
	# file.close()



	####### script for debugging purpose only######
	###############################################
	# start_index = 0
	# end_index = 2
	# print(data_file[0])
	# X_train_batch = []
	# y_train_batch = []
	# # # X_train_batch[start_index:end_index] = np.asarray([create_rxn_Morgan2FP(x[0].split('.'),x[1].split('.')) for x in data_file[start_index:end_index]])
	# # # y_train_batch[start_index:end_index] = np.asarray([x[2] for x in data_file[start_index:end_index]])
	# # # X_train_batch[start_index:end_index] = [a.astype('float32') for a in X_train_batch[start_index:end_index]]
	# # # y_train_batch[start_index:end_index] = [a.astype('float32') for a in y_train_batch[start_index:end_index]]
			
	# # # print(X_train_batch[0])
	# # # print(X_train_batch[0].shape)


	# X_train_batch = []
	# y_train_batch = []
	# for a in batch_data_generator(data_file, 0, 100, batch_size):
	# 	print(type(a))
