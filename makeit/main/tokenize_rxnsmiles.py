from __future__ import print_function
from keras.preprocessing.text import Tokenizer # text pre-processing
import datetime
import cPickle
import json
import sys  # for commanad line
import os

def get_tokenizer_fpath():
	'''Returns file path where tokenizer is backed up'''
	return os.path.join(os.path.dirname(os.path.dirname(__file__)), 'models')

def build_tokenizer(data_fname, N = None):
	'''This function trains a keras.preprocessing.text.Tokenizer on a set
	of data in .json form, where the .json file must be a list of elements 
	where the first value is a reaction smiles string'''

	# Initialize tokenizer
	tokenizer = Tokenizer(nb_words = N, filters = '', lower = True, 
		split = ' ')

	# Load data
	data_fid = open(data_fname, 'r')
	data = json.load(data_fid)
	print('...loaded data from {}'.format(data_fname))
	rxn_smiles_strings = [x[0] for x in data]

	# Iterate through now to build list
	tokenizer.fit_on_texts(rxn_smiles_strings)
	print '...fit tokenizer on reaction smiles strings'

	# Print example
	name = data[0][0]
	print('EXAMPLE:')
	print('    ' + name)
	print('    ' + str(tokenizer.texts_to_sequences([name])))

	# Return trained tokenizer
	return tokenizer

def save_tokenizer(tokenizer, data_fname):
	'''Saves tokenizer object and associated information'''

	# Create filename
	if tokenizer.nb_words:
		fpath_root = os.path.join(get_tokenizer_fpath(), 'tokenizer_rxnsmiles_{}'.format(tokenizer.nb_words))
	else:
		fpath_root = os.path.join(get_tokenizer_fpath(), 'tokenizer_rxnsmiles_{}'.format(len(tokenizer.word_counts)))

	# Dump data
	cPickle.dump(tokenizer, open(fpath_root + '.cpickle', 'wb'))

	# Write to info file
	info_fid = open(fpath_root + '.info', 'w')
	time_now = datetime.datetime.utcnow()
	info_fid.write('{} generated at UTC {}\n\n'.format(fpath_root, time_now))
	info_fid.write('File details\n------------\n')
	info_fid.write('- data source: {}\n'.format(data_fname))
	info_fid.write('- document count: {}\n'.format(tokenizer.document_count))
	info_fid.write('- vocabulary size: {}\n'.format(len(tokenizer.word_counts)))
	info_fid.write('- nb_words: {}\n'.format(tokenizer.nb_words))
	info_fid.close()

	print('...saved tokenizer to {}'.format(fpath_root))
	return True

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print('Usage: {} "data_file.json" [max # vocab]'.format(sys.argv[0]))
		print('    data_file.json must have reaction smiles as first elements')
		quit(1)

	# Build
	if len(sys.argv) == 3:
		tokenizer = build_tokenizer(sys.argv[1], int(sys.argv[2]))
	else:
		tokenizer = build_tokenizer(sys.argv[1])

	# Save
	save_tokenizer(tokenizer, sys.argv[1])