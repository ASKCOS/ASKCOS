# Defines Keras layer class for undirected attributed graph embedding
# based on description in http://arxiv.org/pdf/1509.09292v2.pdfv
'''

Now intended for keras 1.0, not 0.3

'''

import numpy as np
import keras.backend as K
import theano.tensor as T # should write custom back-end eventually, but this is quick fix
import theano
from keras import activations, initializations, regularizers
from keras.engine.topology import Layer
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
		inner_dim: equivalent to F in Duvenaud's paper. the number of attributes
			for each (bond, atom) pair concatenated. Does NOT include the extract
			is_bond_present flag. 
		depth: radius of fingerprint (how many times to recursively mix attributes)
		init_output: initialization for weights in output layer
		activation_output: activation function for output layer. Softmax is recommended
			because it can help increase sparsity, making it more like a real fingerprint
		init_inner: initialization for inner weights for mixing attributes. Identity is
			recommended for the initialization for simplicity
		activation_inner: activation function for the inner layer.
		scale_output: scale to use for output weight initializations. Large output weights
			are closer to a true sparse fingerprint, but small output weights might
			be better to not get stuck in local minima (with low gradients)
		padding: whether to look for padding in the input tensors. 		
	'''

	def __init__(self, output_dim, inner_dim, depth = 2, init_output='uniform', 
			activation_output='softmax', init_inner='identity',
			activation_inner='linear', inner_regularizer=None, 
			output_regularizer=None, scale_output=0.01, padding=False, **kwargs):
		if depth < 1:
			quit('Cannot use GraphFP with depth zero')
		self.init_output = initializations.get(init_output)
		self.activation_output = activations.get(activation_output)
		self.init_inner = initializations.get(init_inner)
		self.activation_inner = activations.get(activation_inner)
		self.output_dim = output_dim
		self.inner_dim = inner_dim
		self.depth = depth
		self.scale_output = scale_output
		self.padding = padding

		# Regularizers (currently broken)
		self.inner_regularizer = regularizers.get(inner_regularizer)
		self.output_regularizer = regularizers.get(output_regularizer)

		self.initial_weights = None
		self.input_dim = 4 # each entry is a 3D N_atom x N_atom x N_feature tensor
		if self.input_dim:
			kwargs['input_shape'] = (None, None, None,) # 3D tensor for each input
		#self.input = K.placeholder(ndim = 4)
		super(GraphFP, self).__init__(**kwargs)


	def build(self, input_shape):
		'''Builds internal weights and paramer attribute'''
		# NOTE: NEED TO TILE AND EVALUATE SO THAT PARAMS CAN BE VARIABLES
		# OTHERWISE K.GET_VALUE() DOES NOT WORK

		# Define template weights for inner FxF
		W_inner = self.init_inner((self.inner_dim, self.inner_dim))
		b_inner = K.zeros((1, self.inner_dim))
		# Initialize weights tensor 
		self.W_inner = K.variable(T.tile(W_inner, (self.depth + 1, 1, 1)).eval() + \
			initializations.uniform((self.depth + 1, self.inner_dim, self.inner_dim)).eval())
		self.W_inner.name = 'T:W_inner'
		self.b_inner = K.variable(T.tile(b_inner, (self.depth + 1, 1, 1)).eval()  + \
			initializations.uniform((self.depth + 1, 1, self.inner_dim)).eval())
		self.b_inner.name = 'T:b_inner'
		# # Concatenate third dimension (depth) so different layers can have 
		# # different weights. Now, self.W_inner[#,:,:] corresponds to the 
		# # weight matrix for layer/depth #.

		# Define template weights for output FxL
		W_output = self.init_output((self.inner_dim, self.output_dim), scale = self.scale_output)
		b_output = K.zeros((1, self.output_dim))
		# Initialize weights tensor
		self.W_output = K.variable(T.tile(W_output, (self.depth + 1, 1, 1)).eval())
		self.W_output.name = 'T:W_output'
		self.b_output = K.variable(T.tile(b_output, (self.depth + 1, 1, 1)).eval())
		self.b_output.name = 'T:b_output'
		# # Concatenate third dimension (depth) so different layers can have 
		# # different weights. Now, self.W_output[#,:,:] corresponds to the 
		# # weight matrix for layer/depth #.

		# # Get regularizers
		# self.regularizers = []
		# if self.inner_regularizer:
		# 	self.inner_regularizer.set_param(self.W_inner)
		# 	self.regularizers.append(self.inner_regularizer)
		# if self.output_regularizer:
		# 	self.output_regularizer.set_param(self.W_output)
		# 	self.regularizers.append(self.output_regularizer)

		# Pack params
		self.trainable_weights = [self.W_inner, 
					   self.b_inner,
					   self.W_output,
					   self.b_output]
		self.params = [self.W_inner, 
					   self.b_inner,
					   self.W_output,
					   self.b_output]

	def get_output_shape_for(self, input_shape):
		return (input_shape[0], self.output_dim)

	def call(self, x, mask=None):
		(output, updates) = theano.scan(lambda x_one: self.get_output_singlesample(x_one), sequences = x)
		return output

	def get_output_singlesample(self, original_graph, train=False):
		'''For a 3D tensor, get the output. Avoids the need for even more complicated vectorization'''
		# Check padding
		if self.padding:
			rowsum = original_graph.sum(axis = 0) # add across
			trim = rowsum[:, -1] # last feature == bond flag
			trim_to = T.eq(trim, 0).nonzero()[0][0] # first index with no bonds
			original_graph = original_graph[:trim_to, :trim_to, :] # reduced graph

		# Get attribute values for layer zero
		# where attributes is a 2D tensor and attributes[#, :] is the vector of
		# concatenated node and edge attributes. In the first layer (depth 0), the 
		# edge attribute section is initialized to zeros. After increaseing depth, howevevr,
		# this part of the vector will become non-zero.

		# The first attributes matrix is just graph_tensor[i, i, :], but we can't use that 
		# kind of advanced indexing
		# Want to extract tensor diagonal as matrix, but can't do that directly...
		# Want to loop over third dimension, so need to dimshuffle
		# print('original graph shape: {}'.format(original_graph.shape.eval()))
		(attributes, updates) = theano.scan(lambda x: x.diagonal(), sequences = original_graph.dimshuffle((2, 0, 1)))
		attributes.name = 'attributes'
		# Now the attributes is (N_features x N_nodes), so we need to transpose
		attributes = attributes.T
		attributes.name = 'attributes post-transpose'

		# Get initial fingerprint
		presum_fp = self.attributes_to_fp_contribution(attributes, 0)
		fp = K.sum(presum_fp, axis = 0) # sum across atom contributions
		fp.name = 'initial fingerprint'

		# Get bond matrix
		bonds = original_graph[:, :, -1] # flag if the bond is present, (N_atom x N_atom)
		bonds.name = 'bonds'

		# Iterate through different depths, updating attributes each time
		for depth in range(self.depth):
			# print('depth {} fp: {}'.format(depth, fp.eval()))
			(attributes, graph) = self.attributes_update(attributes, depth + 1, original_graph, original_graph, bonds)
			presum_fp_new = self.attributes_to_fp_contribution(attributes, depth)
			presum_fp_new.name = 'presum_fp_new contribution'
			fp = fp + K.sum(presum_fp_new, axis = 0) 

		# print('final fp: {}'.format(fp.eval()))
		return fp

	def get_output_singlesample_detail(self, train = False):
		'''For a single input, get the output fingerprint *before* summing by atom and 
		over different depths. This should only be used in visualization.'''
		original_graph = self.get_input(train)[0]
		# Check padding
		if self.padding:
			rowsum = original_graph.sum(axis = 0) # add across
			trim = rowsum[:, -1] # last feature == bond flag
			trim_to = T.eq(trim, 0).nonzero()[0][0] # first index with no bonds
			original_graph = original_graph[:trim_to, :trim_to, :] # reduced graph
		(attributes, updates) = theano.scan(lambda x: x.diagonal(), sequences = original_graph.dimshuffle((2, 0, 1)))
		# Now the attributes is (N_features x N_nodes), so we need to transpose
		attributes = attributes.T
		# Get initial fingerprint (not summed over atoms)
		fp = self.attributes_to_fp_contribution(attributes, 0)
		# Pad right (since this will become tensor)
		fp = T.shape_padright(fp)
		# Get bond matrix
		bonds = original_graph[:, :, -1] # flag if the bond is present, (N_atom x N_atom)
		# Iterate through different depths, updating attributes each time
		for depth in range(self.depth):
			depth = depth + 1 # correct for zero-indexing
			(attributes, graph) = self.attributes_update(attributes, depth, original_graph, original_graph, bonds)
			fp_new = self.attributes_to_fp_contribution(attributes, depth)
			fp = K.concatenate((fp, T.shape_padright(fp_new)), axis = 2)
		return fp


	def attributes_update(self, attributes, depth, graph, original_graph, bonds):
		'''Given the current attributes, the current depth, and the graph that the attributes
		are based on, this function will update the 2D attributes tensor'''
		
		############# GET NEW ATTRIBUTE MATRIX #########################
		# New pre-activated attribute matrix v = M_i,j,: x ones((N_atom, 1)) -> (N_atom, N_features) 
		# as long as dimensions are appropriately shuffled
		shuffled_graph = graph.copy().dimshuffle((2, 0, 1)) # (N_feature x N_atom x N_atom)
		shuffled_graph.name = 'shuffled_graph'
		# print('shuffled shape of graph: {}'.format(shuffled_graph.eval().shape))
		ones_vec = K.ones_like(attributes[:, 0]) # (N_atom x 1)
		ones_vec.name = 'ones_vec'
		

		# Embed individually
		# (scan sequences iterates over the FIRST dimension)
		# (flatten(ndim) keeps the first ndim-1 dimensions the same, then expands the rest to fill)
		flattened_graph = shuffled_graph.flatten(ndim = 2).T # (N_atom^2 x N_feature)
		# Embed each possible atom-atom interaction
		(new_presummed_attributes_flat, updates) = theano.scan(lambda x: self.activation_inner(
				K.dot(x[:-1], self.W_inner[depth, :, :]) + self.b_inner[depth, 0, :]), sequences = flattened_graph) # still (N_atom^2 x N_feature)
		# Reshape into #(N_feature-1 x N_atom x N_atom)
		new_presummed_attributes = new_presummed_attributes_flat.T.reshape(shuffled_graph[:-1,:,:].shape)

		# Now sum activated self+neighbors
		(new_attributes, updates) = theano.scan(lambda x: K.dot(x, ones_vec), sequences = new_presummed_attributes) # (N_features x N_atom)
		# print('SIZE OF PRE_ACTIVATED_ATTRIBUTES: {}'.format(new_preactivated_attributes.eval().shape))

		# Append last feature (bond flag) after the loop
		# print('SIZE OF NEW_ATTRIBUTES BEFORE CONCAT: {}'.format(new_attributes.eval().shape))
		new_attributes = K.concatenate((new_attributes.T, attributes[:, -1:]), axis = 1)
		new_attributes.name = 'new_attributes'
		# print('new_attributes shape: {}'.format(new_attributes.eval().shape))


		############ UPDATE GRAPH TENSOR WITH NEW ATOM ATTRIBUTES ###################
		### Node attribute contribution is located in every entry of graph[i,j,:] where
		### there is a bond @ ij or when i = j (self)
		# Get atoms matrix (identity)
		atoms = T.identity_like(bonds) # (N_atom x N_atom)
		atoms.name = 'atoms_identity'
		# Combine
		bonds_or_atoms = bonds + atoms # (N_atom x N_atom)
		bonds_or_atoms.name = 'bonds_or_atoms'
		# print('shape of bonds_or_atoms: {}'.format(bonds_or_atoms.eval().shape))

		atom_indeces = T.arange(ones_vec.shape[0]) # 0 to N_atoms - 1 (indeces)
		atom_indeces.name = 'atom_indeces vector'
		### Subtract previous node attribute contribution
		# Multiply each entry in bonds_or_atoms by the previous atom features for that column
		(old_features_to_sub, updates) = theano.scan(lambda i: T.outer(bonds_or_atoms[:, i], attributes[i, :]), 
			sequences = T.arange(ones_vec.shape[0]))
		old_features_to_sub.name = 'old_features_to_sub'
		# print('old_to_sub: {}'.format(old_features_to_sub.eval()))

		### Add new node ttribute contribution
		# Multiply each entry in bonds_or_atoms by the previous atom features for that column
		(new_features_to_add, updates) = theano.scan(lambda i: T.outer(bonds_or_atoms[:, i], new_attributes[i, :]),
			sequences = T.arange(ones_vec.shape[0]))
		new_features_to_add.name = 'new_features_to_add'
		# print('new_to_add: {}'.format(new_features_to_add.eval()))

		# Update new graph
		new_graph = graph - old_features_to_sub + new_features_to_add
		new_graph.name = 'new_graph'

		return (new_attributes, new_graph)


	def attributes_to_fp_contribution(self, attributes, depth):
		'''Given a 2D tensor of attributes where the first dimension corresponds to a single
		node, this method will apply the output sparsifying (often softmax) function and return
		the contribution to the fingerprint'''
		# Apply output activation function
		# print('SIZE OF ATTRIBUTES[:, :-1]: {}'.format(attributes[:, :-1].eval().shape))
		output_dot = K.dot(attributes[:, :-1], self.W_output[depth, :, :]) # ignore last attribute (bond flag)
		output_dot.name = 'output_dot'
		# print('OUTPUT DOT')
		# print(output_dot.eval())
		output_bias = self.b_output[depth, 0, :]
		output_bias.name = 'output_bias'
		# print('OUTPUT BIAS')
		# print(output_bias.eval())
		output_activated = self.activation_output(output_dot + output_bias)
		output_activated.name = 'output_activated'
		return output_activated

	def get_config(self):
		config = {'output_dim': self.output_dim,
				  'inner_dim' : self.inner_dim,
				  'init_output' : self.init_output.__name__,
				  'init_inner' : self.init_inner.__name__,
				  'activation_inner': self.activation_inner.__name__,
				  'activation_output' : self.activation_output.__name__,
				  'inner_regularizer' : self.inner_regularizer.get_config() if self.inner_regularizer else None,
				  'output_regularizer' : self.output_regularizer.get_config() if self.output_regularizer else None,
				  'input_dim': self.input_dim,
				  'depth' : self.depth}
		base_config = super(GraphFP, self).get_config()
		return dict(list(base_config.items()) + list(config.items()))
