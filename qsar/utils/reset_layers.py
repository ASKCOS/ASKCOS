from __future__ import print_function

from keras import initializations
import numpy as np
import theano.tensor as T 

def reset(model):
	'''Given a Keras model consisting only of GraphFP, Dense, and Dropout layers,
	this function will reset the trainable weights to save time for CV tests.'''

	for layer in model.layers:
		# Note: these are custom depending on the layer type
		if '.GraphFP' in str(layer):
			W_inner = layer.init_inner((layer.inner_dim, layer.inner_dim))
			b_inner = np.zeros((1, layer.inner_dim))
			# Inner weights
			layer.W_inner.set_value((T.tile(W_inner, (layer.depth + 1, 1, 1)).eval() + \
				initializations.uniform((layer.depth + 1, layer.inner_dim, layer.inner_dim)).eval()).astype(np.float32))
			layer.b_inner.set_value((T.tile(b_inner, (layer.depth + 1, 1, 1)).eval()  + \
				initializations.uniform((layer.depth + 1, 1, layer.inner_dim)).eval()).astype(np.float32))

			# Outer weights
			W_output = layer.init_output((layer.inner_dim, layer.output_dim), scale = layer.scale_output)
			b_output = np.zeros((1, layer.output_dim))
			# Initialize weights tensor
			layer.W_output.set_value((T.tile(W_output, (layer.depth + 1, 1, 1)).eval()).astype(np.float32))
			layer.b_output.set_value((T.tile(b_output, (layer.depth + 1, 1, 1)).eval()).astype(np.float32))
			print('graphFP layer reset')

		elif '.Dense' in str(layer):
			layer.W.set_value((layer.init(layer.W.shape.eval()).eval()).astype(np.float32))
			layer.b.set_value(np.zeros(layer.b.shape.eval(), dtype=np.float32))
			print('dense layer reset')

		elif '.Dropout' in str(layer):
			print('dropout unchanged')
		else:
			raise ValueError('Unknown layer {}, cannot reset weights'.format(str(layer)))
	print('Reset model weights')
	return model