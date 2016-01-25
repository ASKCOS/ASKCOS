from __future__ import print_function
from keras.preprocessing.text import Tokenizer # text pre-processing
import collections
import datetime
import cPickle
import json
import sys  # for commanad line
import os

def get_tokenizer_fpath():
	'''Returns file path where tokenizer is backed up'''
	return os.path.join(os.path.dirname(os.path.dirname(__file__)), 'models')

def build_tokenizer(data_fname, N = None):
	'''This function defines a dictionary character mapping'''

	# Load data
	data_fid = open(data_fname, 'r')
	data = json.load(data_fid)
	print('...loaded data from {}'.format(data_fname))
	rxn_smiles_strings = [str(x[0].strip()) for x in data]

	# Merge
	all_chars = ''.join(rxn_smiles_strings)
	print('total characters read: {}'.format(len(all_chars)))

	# Define collection
	chars = collections.Counter(all_chars)

	# Get unique values ordered by frequency
	chars_by_freq = sorted(chars, key = chars.get, reverse = True)

	tokenizer = {}
	for (i, char) in enumerate(chars_by_freq):
		tokenizer[char] = i + 1 # save 0 for masking

	# Print example
	name = rxn_smiles_strings[0]
	print('EXAMPLE:')
	print('    {}'.format(name))
	print('    {}'.format([tokenizer[x] for x in name]))

	# Return trained tokenizer
	return tokenizer

def save_tokenizer(tokenizer, data_fname):
	'''Saves tokenizer object and associated information'''

	# Create filename
	fpath_root = os.path.join(get_tokenizer_fpath(), 'tokenizer_rxnsmiles')

	# Dump data
	json.dump(tokenizer, open(fpath_root + '.json', 'wb'))

	# Write to info file
	info_fid = open(fpath_root + '.info', 'w')
	time_now = datetime.datetime.utcnow()
	info_fid.write('{} generated at UTC {}\n\n'.format(fpath_root, time_now))
	info_fid.write('File details\n------------\n')
	info_fid.write('- data source: {}\n'.format(data_fname))
	info_fid.write('- vocabulary size: {}\n'.format(len(tokenizer.keys())))
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