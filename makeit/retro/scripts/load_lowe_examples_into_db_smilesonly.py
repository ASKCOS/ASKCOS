## Meant to use with Daniel Lowe's patent database
# reaction SMILES strings only (for now)

from __future__ import print_function
import argparse
from numpy.random import shuffle # for random selection
from numpy.random import random
import rdkit.Chem as Chem          # molecule building
from rdkit.Chem import AllChem
from collections import defaultdict
import rdkit.Chem.Draw as Draw
from rdkit import RDLogger
import datetime # for info files
import json # for dumping
import sys  # for commanad line
import os   # for file paths
import re 
import itertools

# DATABASE
from pymongo import MongoClient    # mongodb plugin
client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaction_examples']
collection = db['lowe_1976-2013_USPTOgrants']

def main(db_fpath, N = 15):
	'''Read reactions from Lowe's patent reaction SMILES'''

	# Open file
	data_fid = open(db_fpath, 'r')

	# Define scoring variables
	total_templates = 0 # total reactions simulated (excludes skipped)
	total_correct = 0 # actual products predicted
	total_precise = 0 # ONLY actual products predicted

	try: # to allow breaking
		# Look for entries
		documents = []

		for i, line in enumerate(data_fid):
			# Are we done?
			if i == N:
				break

			# Unpack
			line = line.strip()
			reaction_smiles = line.split('\t')[0]
			reference = line.split('\t')[1]
			reaction_smiles = reaction_smiles.split(' ')[0]

			# Load into database
			documents.append(
				{
					'reaction_smiles': reaction_smiles,
					'reference': reference,
					'random': random(),
				}
			)

			# Report progress and insert every 1000
			if ((i+1) % 1000) == 0:
				print('{}/{}'.format(i+1, N))
				result = collection.insert(documents)
				documents = []
		result = collection.insert(documents)
	except KeyboardInterrupt:
		print('Stopped early!')		
	except Exception as e:
		print(e)
		print(line)

	print('Created {} database entries'.format(collection.find().count()))

	return True


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('data_file', type = str, 
		 				help = 'File where each line is an atom-mapped smiles reaction')
	parser.add_argument('-n', '--num', type = int, default = 50,
						help = 'Maximum number of records to load; defaults to 50')
	args = parser.parse_args()

	clear = raw_input('Do you want to clear the {} existing examples? '.format(collection.find().count()))
	if clear in ['y', 'Y', 'yes', '1', 'Yes']:
		result = collection.delete_many({})
		print('Cleared {} entries from collection'.format(result.deleted_count))

	main(args.data_file, N = args.num)