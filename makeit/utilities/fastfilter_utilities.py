
#import files
from __future__ import division

import theano
from theano.tensor import lt,le,eq,gt,ge

import numpy as np#something

from keras import backend as K 
from keras.layers import Dense, Activation, Input
from keras.layers import merge, activations
from keras.optimizers import SGD, Adam, Adadelta
#import keras.engine.topology 
from keras.engine.topology import Layer
import os


# lg = RDLogger.logger()
# lg.setLevel(RDLogger.CRITICAL)
####utilities
def set_keras_backend(backend):

    if K.backend() != backend:
        os.environ['KERAS_BACKEND'] = backend
        reload(K)
        assert K.backend() == backend

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


def pos_ct(y_true, y_pred):
	pos_pred = K.sum(gt((K.clip(y_pred, 0, 1)),0.5))
	return pos_pred
def true_pos(y_true, y_pred):
	true_pos_ct = K.sum(gt((K.clip(y_pred*y_true, 0, 1)),0.5))
	return true_pos_ct
def real_pos(y_true, y_pred):
	real_pos_ct = K.sum(gt((K.clip(y_true, 0, 1)),0.5))
	return real_pos_ct

