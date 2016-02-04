# Defines Keras layer class for undirected attributed graph embedding
# based on description in http://arxiv.org/pdf/1509.09292v2.pdfv

import numpy as np
import keras.backend as K
from keras import activations, initializations
from keras.layers.core import Layer

class GraphFP(Layer):
	'''Embedding layer for undirected, attributed graphs following the ideas of the 
	extended connectivity fingerprints (ECFPs) or functional connectivity FPs (FCFPs).
	It must be the first layer in a model.

	# Input shape
		1D array with shape: `(nb_samples, 1)`.

	# Output shape
		2D tensor with shape: `(nbb_samples, output_dim)`.

	# Arguments
		output_dim: int > 0, size of the fingerprint
		init: name of initialization function for the weights of the layer
			(see [initializations](keras/initializations.md)),
			or alternatively, Theano function to use for weights
			initialization. This parameter is only relevant
			if you don't pass a `weights` argument.
		activation: activation function used for the *inner* smoothing function,
			the outer layer always uses softmax.
		inner_dim: equivalent to F in Duvenaud's paper. the number of attributes
			for each (bond, atom) pair concatenated.
	'''

	def __init__(self, output_dim, inner_dim, depth = 2, init_output='glorot_uniform', 
			activation_output='softmax', init_inner='identity',
			activation_inner='linear', activity_regularizer=None, **kwargs):
		self.init_output = initializations.get(init_output)
		self.activation_output = activations.get(activation_output)
		self.init_inner = initializations.get(init_inner)
		self.activation_inner = activations.get(activation_inner)
		self.output_dim = output_dim
		self.inner_dim = inner_dim
		self.depth = depth

		self.initial_weights = None
		self.input_dim = 1
		if self.input_dim:
			kwargs['input_shape'] = (self.input_dim,)
		self.input = K.placeholder(ndim=1)
		super(GraphFP, self).__init__(**kwargs)

	def build(self):

		# Define template weights for inner FxF
		W_inner = self.init_inner((self.inner_dim, self.inner_dim))
		b_inner = K.zeros((1, self.inner_dim))
		# Initialize weights tensor 
		self.W_inner = W_inner.copy()
		self.b_inner = b_inner.copy()
		# Concatenate third dimension (depth) so different layers can have 
		# different weights. Now, self.W_inner[:,:,#] corresponds to the 
		# weight matrix for layer/depth #.
		for i in range(self.depth):
			self.W_inner = np.dstack((self.W_inner, W_inner.copy()))
			self.b_inner = np.dstack((self.b_inner, b_inner.copy()))

		# Define template weights for output FxL
		W_output = self.init_output((self.inner_dim, self.output_dim))
		b_output = K.zeros((1, self.output_dim))
		# Initialize weights tensor
		self.W_output = W_output.copy()
		self.b_output = b_output.copy()
		# Concatenate third dimension (depth) so different layers can have 
		# different weights. Now, self.W_output[:,:,#] corresponds to the 
		# weight matrix for layer/depth #.
		for i in range(self.depth):
			self.W_output = np.dstack((self.W_output, W_output.copy()))
			self.b_output = np.dstack((self.b_output, b_output.copy()))

		self.params = [self.W_inner, 
					   self.b_inner,
					   self.W_output,
					   self.b_output]

	@property
	def output_shape(self):
		return (self.input_shape[0], self.output_dim)

	def get_output(self, train=False):
		graph = self.get_input(train)

		# Get attribute values for layer zero
		# where attributes is a 2D tensor and attributes[#, :] is the vector of
		# concatenated node and edge attributes. In the first layer (depth 0), the 
		# edge attribute section is initialized to zeros. After increaseing depth, howevevr,
		# this part of the vector will become non-zero.
		attributes = K.concatenate((graph.nodeAttributes(), K.zeros_like(graph.edgeAttributes())), axis = 1)
		fp = self.attributes_to_fp_contribution(attributes, 0)

		# Iterate through different depths, updating attributes each time
		for depth in range(self.depth):
			depth = depth + 1 # correct for zero-indexing

			attributes = self.attributes_update(attributes, depth, graph)
			fp += self.attributes_to_fp_contribution(attributes, depth)

		return fp

	def attributes_update(self, attributes, depth, graph):
		'''Given the current attributes, the current depth, and the graph that the attributes
		are based on, this function will update the 2D attributes tensor'''
		# Get original attribute vectors for sizing
		edgeAttributes = graph.edgeAttributes()
		nodeAttributes = graph.nodeAttributes()
		# Copy new attributes so we can overwrite
		new_attributes = attributes.copy()
		for i, node in enumerate(graph.nodes):
			v = attributes[i, :].copy() # initialize with current attributes
			for (ni, ei) in node.neighbors:
				# mix in contributions from neighbor node's attributes
				v += attributes[ni, :] 
				# modify with edge information to get to that neighbor
				v += K.concatenate((K.zeros_like(nodeAttributes[0]), edgeAttributes[ei]), axis = 1)
			# smooth with inner activation function and weights
			new_attributes[i, :] = self.activation_inner(
				K.dot(v, self.W_inner[:, :, depth]) + 
				self.b_inner
			)
  		# Return updated attribute tensor
		return new_attributes


	def attributes_to_fp_contribution(self, attributes, depth):
		'''Given a 2D tensor of attributes where the first dimension corresponds to a single
		node, this method will apply the output sparsifying (often softmax) function and return
		the contribution to the fingerprint'''
		# Apply output activation function
		return K.sum(
					self.activation_output(
						K.dot(attributes, self.W_output[:, :, depth]) + 
						K.ones((K.shape(attributes)[0], 1)) * self.b_output[:, :, depth]
					)
				)

	def get_config(self):
		config = {'name': self.__class__.__name__,
				  'output_dim': self.output_dim,
				  'inner_dim' : self.inner_dim,
				  'init_output' : self.init_output.__name__,
				  'init_inner' : self.init_inner.__name__,
				  'activation_inner': self.activation_inner.__name__,
				  'activation_output' : self.activation_output.__name__,
				  'input_dim': self.input_dim,
				  'depth' : self.depth}
		base_config = super(GraphFP, self).get_config()
		return dict(list(base_config.items()) + list(config.items()))
