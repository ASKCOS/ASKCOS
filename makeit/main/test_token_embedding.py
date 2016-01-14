from makeit.utils.chemnet_connect import * # mongodb connection, gets 'chemicals', 'reactions'
from makeit.main.tokenize_chemicals import load_tokenizer
from makeit.utils.parsing import SplitChemicalName 
from makeit.utils.database import Collection_Manager
from keras.models import Sequential
from keras.layers.embeddings import Embedding
from keras.optimizers import SGD
import cPickle
import os

model_fname = os.path.join(os.path.dirname(__file__), 'test_model.cpickle')

def build_model(input_dim, output_dim):
	'''Generates simple embedding model to use tokenized chemical name as
	input in order to predict a single-valued output (i.e., MW)'''

	# Base model
	model = Sequential()

	# Add embedding layer
	model.add(Embedding(input_dim, output_dim, 
		init = 'uniform', input_length = None, W_regularizer = None, 
		activity_regularizer = None, W_constraint = None, 
		mask_zero = False, weights = None))
	
	# Define optimizer
	sgd = SGD(lr = 0.1, decay = 1e-6, momentum = 0.9)

	# Compile
	model.compile(loss = 'mean_squared_error', optimizer = sgd)

	return model

def save_model(model):
	'''Saves NN model object according to the filename defined
	in makeit.main.tokenize_chemicals.py'''
	cPickle.dump(model, open(model_fname, 'wb'))
	return True

def load_model():
	'''Loads NN model object according to the filename defined
	in makeit.main.tokenize_chemicals.py'''
	return cPickle.load(open(model_fname, 'rb'))

def train_model(model, tokenizer, chemicals):
	'''Trains the model to predict chemical molecular weights'''

	# Put wrapper around chemicals collection
	Chemicals = Collection_Manager(chemicals, batch_size = 100)

	# Look through database
	#while not Chemicals.done:
	for i in range(10):
		# Pull a batch of documents for training
		batch_docs  = Chemicals.pull()
		batch_names = []
		batch_mws = []
		for doc in batch_docs:
			if doc['name'] and doc['mol_weight']:
				batch_names.append(' '.join(SplitChemicalName(doc['name'][0])))
				batch_mws.append(doc['mol_weight'])

		# Run through tokenizer
		batch_matrix = tokenizer.texts_to_matrix(batch_names)

		# Train on the batch
		model.train_on_batch(batch_matrix, batch_mws)
		
	# Evaluate using first batch?
	Chemicals.index = 0
	batch_docs  = Chemicals.pull()
	batch_names = [x['name'][0] for x in batch_docs]
	batch_names = [' '.join(SplitChemicalName(x)) for x in batch_names]
	batch_mws   = [x['mol_weight'] for x in batch_docs]
	batch_matrix = tokenizer.texts_to_matrix(batch_names)
	score = model.evaluate(batch_matrix, batch_mws)

	return (model, score)


if __name__ == '__main__':

	# Did we already fit tokenizer?
	try:
		# Load
		tokenizer = load_tokenizer()
	# Need to train it
	except:
		# Don't use this module to do that
		raise IOError('Could not load tokenizer.cpickle, ' +
				'try running makeit.main.tokenize_chemicals.py ' +
				'to fit tokenizer to chemical vocabulary')

	# Did we already define model?
	try:
		# Load
		model = load_model()
		print 'restored model from %s' % model_fname
	# Need to build it
	except:
		# Build model
		model = build_model(len(tokenizer.word_counts), 1)
		print 'built untrained model'
		# Save for future
		save_model(model)
		print 'saved model'
		
	# Train model
	(model, score) = train_model(model, tokenizer, chemicals)
	print 'trained model'
	print '-> score = ' + str(score)
		

