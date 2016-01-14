from makeit.utils.chemnet_connect import * # mongodb connection, gets 'chemicals', 'reactions'
from makeit.utils.parsing import SplitChemicalName 
from keras.preprocessing.text import Tokenizer # text pre-processing
from random import randint
import cPickle
import os

tokenizer_fname = os.path.join(os.path.dirname(__file__), 'tokenizer.cpickle')

def build_tokenizer(N = 10000):
	'''This function trains a keras.preprocessing.text.Tokenizer on N
	randomly chosen chemicals from the chemnet database *with* replacement,
	which should be fine for low 'N' '''
	# Initialize tokenizer
	tokenizer = Tokenizer(nb_words = None, filters = '', lower = True, 
		split = ' ')

	# IGeneratae random index numbers for database
	Ntot = chemicals.find().count()
	db_indeces = [randint(0, N - 1) for x in range(N)]

	# Iterate through now and pull records
	for i in range(N):

		# Pull record randomly (with replacement)
		chemical = chemicals.find()[db_indeces[i]]

		# Train tokenizer
		if chemical['name']:
			for name in chemical['name']:
				tokenizer.fit_on_texts(SplitChemicalName(name))

		if (i % 1000) == 0:
			print '...reading %i/%i' % (i, N)

	# Return trained tokenizer
	return tokenizer

def save_tokenizer(tokenizer):
	'''Saves tokenizer object according to the filename defined
	in makeit.main.tokenize_chemicals.py'''
	cPickle.dump(tokenizer, open(tokenizer_fname, 'wb'))
	return True

def load_tokenizer():
	'''Loads tokenizer object according to the filename defined
	in makeit.main.tokenize_chemicals.py'''
	return cPickle.load(open(tokenizer_fname, 'rb'))

if __name__ == '__main__':

	# Did we already fit tokenizer?
	try:
		# Load
		tokenizer = load_tokenizer()
		print 'restored tokenizer from %s' % tokenizer_fname
	# Need to train it
	except:
		# Build tokenizer
		tokenizer = build_tokenizer(N = 10000)
		# Save for future
		save_tokenizer(tokenizer)

	# Test with a simple case
	test_doc = chemicals.find_one()
	name = test_doc['name'][0]
	name = ' '.join(SplitChemicalName(name))
	print 'EXAMPLE:'
	print '    ' + name
	print '    ' + str(tokenizer.texts_to_sequences([name]))
