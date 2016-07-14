# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import numpy as np
from scipy.sparse import coo_matrix
import cPickle as pickle
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import os
import sys
from makeit.embedding.descriptors import rxn_level_descriptors
import time

FROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data_edits')

def get_candidates(n = 2, seed = None, outfile = '.', shuffle = False, skip = 0, padUpTo = 500):
	'''
	Pull n example reactions, their candidates, and the true answer
	'''

	from pymongo import MongoClient    # mongodb plugin
	from rdkit import RDLogger
	lg = RDLogger.logger()
	lg.setLevel(4)
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['prediction']
	examples = db['candidate_edits']

	# Define generator
	class Randomizer():
		def __init__(self, seed):
			self.done_ids = []
			np.random.seed(seed)
			if outfile:
				with open(os.path.join(outfile, 'preprocess_candidate_edits_seed.txt'), 'w') as fid:
					fid.write('{}'.format(seed))
		def get_rand(self):
			'''Random WITHOUT replacement'''
			while True:
				try:
					doc = examples.find({'found': True, \
						'random': { '$gte': np.random.random()}}).sort('random', 1).limit(1)
					if not doc: continue
					if doc[0]['_id'] in self.done_ids: continue
					self.done_ids.append(doc[0]['_id'])
					yield doc[0]
				except KeyboardInterrupt:
					print('Terminated early')
					quit(1)
				except:
					pass

	if seed == None:
		seed = np.random.randint(10000)
	else:
		seed = int(seed)
	randomizer = Randomizer(seed)
	if shuffle:
		generator = enumerate(randomizer.get_rand())
	else:
		generator = enumerate(examples.find({'found': True}))

	# Initialize (this is not the best way to do this...)
	reaction_candidate_edits = []
	reaction_candidate_smiles = []
	reaction_true_onehot = []
	reaction_true = []
	for i, reaction in generator:
		if i < skip: continue
		if i == skip + n: break

		candidate_smiles = [a for (a, b) in reaction['edit_candidates']]
		candidate_edits =    [b for (a, b) in reaction['edit_candidates']]
		reactant_smiles = reaction['reactant_smiles']
		product_smiles_true = reaction['product_smiles_true']

		bools = [product_smiles_true == x for x in candidate_smiles]
		print('rxn. {} : {} true entries out of {}'.format(i, sum(bools), len(bools)))
		if sum(bools) > 1:
			pass
			# print('More than one true?')
			# print(reactant_smiles)
			# for (edit, prod) in [(edit, prod) for (boolean, edit, prod) in zip(bools, candidate_edits, candidate_smiles) if boolean]:
			# 	print(prod)
			# 	print(edit)
			# raw_input('Pausing...')
			# continue
		if sum(bools) == 0:
			print('##### True product not found / filtered out #####')
			continue

		# Sort together and append
		zipsort = sorted(zip(bools, candidate_smiles, candidate_edits))
		zipsort = [[(y, z, x) for (y, z, x) in zipsort if y == 1][0]] + \
				  [(y, z, x) for (y, z, x) in zipsort if y == 0]
		zipsort = zipsort[:padUpTo]
		reaction_candidate_edits.append([x for (y, z, x) in zipsort])
		reaction_true_onehot.append([y for (y, z, x) in zipsort])
		reaction_candidate_smiles.append([z for (y, z, x) in zipsort])
		reaction_true.append(str(reactant_smiles) + '>>' + str(product_smiles_true))

	return reaction_candidate_edits, reaction_true_onehot, reaction_candidate_smiles, reaction_true

if __name__ == '__main__':
	n = 5
	padUpTo = 500
	shuffle = False
	skip = 0
	if len(sys.argv) >= 2:
		n = int(sys.argv[1])
	if len(sys.argv) >= 3:
		skip = int(sys.argv[2])

	reaction_candidate_edits, reaction_true_onehot, reaction_candidate_smiles, reaction_true = \
			get_candidates(n = n, shuffle = shuffle, skip = skip, padUpTo = padUpTo)

	rounded_time = '{}-{}'.format(skip, skip + n - 1)

	with open(os.path.join(FROOT, '{}_candidate_edits.pickle'.format(rounded_time)), 'wb') as outfile:
		pickle.dump(reaction_candidate_edits, outfile, pickle.HIGHEST_PROTOCOL)
	with open(os.path.join(FROOT, '{}_candidate_bools.pickle'.format(rounded_time)), 'wb') as outfile:
		pickle.dump(reaction_true_onehot, outfile, pickle.HIGHEST_PROTOCOL)
	with open(os.path.join(FROOT, '{}_candidate_smiles.pickle'.format(rounded_time)), 'wb') as outfile:
		pickle.dump(reaction_candidate_smiles, outfile, pickle.HIGHEST_PROTOCOL)
	with open(os.path.join(FROOT, '{}_reaction_string.pickle'.format(rounded_time)), 'wb') as outfile:
		pickle.dump(reaction_true, outfile, pickle.HIGHEST_PROTOCOL)