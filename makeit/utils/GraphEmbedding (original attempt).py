# Defines Keras layer class for undirected attributed graph embedding
# based on description in http://arxiv.org/pdf/1509.09292v2.pdfv

import numpy as np
import keras.backend as K
import theano.tensor as T # should write custom back-end eventually, but this is quick fix
import theano
from keras import activations, initializations
from keras.layers.core import Layer
from makeit.utils.neural_fp import Graph, molToGraph
from rdkit.Chem import MolFromSmiles

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
		if depth < 1:
			quit('Cannot use GraphFP with depth zero')
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
		#self.input = K.placeholder(ndim=1)
		super(GraphFP, self).__init__(**kwargs)
		self.input = theano.gof.graph.Variable(molToGraph(MolFromSmiles('CC'))) # override?
		self.input = molToGraph(MolFromSmiles('CC')) # override?
		self.input.name = 'placeholder_input'


	def build(self):
		'''Builds internal weights and paramer attribute'''
		# NOTE: NEED TO TILE AND EVALUATE SO THAT PARAMS CAN BE VARIABLES
		# OTHERWISE K.GET_VALUE() DOES NOT WORK

		# Define template weights for inner FxF
		W_inner = self.init_inner((self.inner_dim, self.inner_dim))
		b_inner = K.zeros((1, self.inner_dim))
		# Initialize weights tensor 
		self.W_inner = K.variable(T.tile(W_inner, (self.depth + 1, 1, 1)).eval())
		self.W_inner.name = 'T:W_inner'
		self.b_inner = K.variable(T.tile(b_inner, (self.depth + 1, 1, 1)).eval())
		self.b_inner.name = 'T:b_inner'
		# # Concatenate third dimension (depth) so different layers can have 
		# # different weights. Now, self.W_inner[#,:,:] corresponds to the 
		# # weight matrix for layer/depth #.

		# Define template weights for output FxL
		W_output = self.init_output((self.inner_dim, self.output_dim))
		b_output = K.zeros((1, self.output_dim))
		# Initialize weights tensor
		self.W_output = K.variable(T.tile(W_output, (self.depth + 1, 1, 1)).eval())
		self.W_output.name = 'T:W_output'
		self.b_output = K.variable(T.tile(b_output, (self.depth + 1, 1, 1)).eval())
		self.b_output.name = 'T:b_output'
		# # Concatenate third dimension (depth) so different layers can have 
		# # different weights. Now, self.W_output[#,:,:] corresponds to the 
		# # weight matrix for layer/depth #.

		# Pack params
		self.params = [self.W_inner, 
					   self.b_inner,
					   self.W_output,
					   self.b_output]

	@property
	def output_shape(self):
		return (self.input_shape[0], self.output_dim)

	def get_output(self, train=False):
		graph = self.get_input(train)
		#graph = graph.eval()
		print(graph)

		# Get attribute values for layer zero
		# where attributes is a 2D tensor and attributes[#, :] is the vector of
		# concatenated node and edge attributes. In the first layer (depth 0), the 
		# edge attribute section is initialized to zeros. After increaseing depth, howevevr,
		# this part of the vector will become non-zero.
		empty_edge_attribute_vec = K.zeros_like(graph.edgeAttributes()[0:1, :])
		empty_edge_attribute_vec.namae = 'empty_edge_attribute_vec'
		empty_edge_attribute_mat = T.tile(empty_edge_attribute_vec, (graph.nodeAttributes().shape[0], 1))
		empty_edge_attribute_mat.name = 'empty_edge_attribute_mat'
		attributes = K.concatenate((graph.nodeAttributes(), empty_edge_attribute_mat), axis = 1)

		# Get initial fingerprint
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
		vs = []
		for i, node in enumerate(graph.nodes):
			v = attributes[i, :].copy() # initialize with current attributes
			v.name = 'v'
			for (ni, ei) in node.neighbors:
				# mix in contributions from neighbor node's attributes
				v += attributes[ni, :] 
				# modify with edge information to get to that neighbor
				v += K.concatenate((K.zeros_like(nodeAttributes[0]), edgeAttributes[ei]), axis = 0)
			# smooth with inner activation function and weights
			inner_dot = K.dot(v, self.W_inner[depth, :, :])
			inner_dot.name = 'inner_dot'
			inner_bias = self.b_inner[depth, 0, :]
			inner_bias.name = 'inner_bias'
			inner_activated = self.activation_inner(inner_dot + inner_bias)
			inner_activated.name = 'inner_activated'
			vs.append(inner_activated)
		# Return updated attribute tensor
		stacked = T.stack(vs)
		stacked.name = 'stacked_attributes'
		return T.stack(vs)


	def attributes_to_fp_contribution(self, attributes, depth):
		'''Given a 2D tensor of attributes where the first dimension corresponds to a single
		node, this method will apply the output sparsifying (often softmax) function and return
		the contribution to the fingerprint'''
		# Apply output activation function
		output_dot = K.dot(attributes, self.W_output[depth, :, :])
		output_dot.name = 'output_dot'
		output_bias = K.ones_like(attributes[:, 0]) * self.b_output[depth, 0, :]
		output_bias.name = 'output_bias'
		output_activated = self.activation_output(output_dot + output_bias)
		output_activated.name = 'output_activated'
		summed_activated = K.sum(output_activated, axis = 0)
		summed_activated.name = 'summed_activated'
		return summed_activated

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
