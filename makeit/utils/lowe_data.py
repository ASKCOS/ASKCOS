# Utilities for using the Daniel Lowe patent data

from __future__ import print_function
import os
from xml.dom import minidom

def get_reaction_file(fpath):
	'''This function is used to recursively traverse a directory and find
	all .cml files contained within. It is used as a generator.'''
	for dir, subdirs, files in os.walk(fpath):
		subdirs.sort(reverse = True)
		for subdir in subdirs:
			for p in get_reaction_file(subdir):
				yield p
		files.sort(reverse = True)
		for file in files:
			yield os.path.join(dir, file)

if __name__ == '__main__':
	root = '/afs/csail.mit.edu/u/c/ccoley/grants'
	print('is_directory? {}'.format(os.path.isdir(root)))
	file_generator = get_reaction_file(root)
	print(file_generator)
	for i, rxn in enumerate(file_generator):
		print('{}: {}'.format(i, rxn))
		document = minidom.parse(rxn)
		rxn_smiles = document.getElementsByTagName('dl:reactionSmiles')[0]
		print(rxn_smiles.firstChild.nodeValue)
		if i == 1:
			break
