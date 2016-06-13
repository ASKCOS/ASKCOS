# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import numpy as np
from scipy.sparse import coo_matrix
import cPickle as pickle
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import os
from makeit.embedding.descriptors import rxn_level_descriptors
from keras.preprocessing.sequence import pad_sequences

def get_candidates(n = 2, seed = None, outfile = '.'):
	'''
	Pull n example reactions, their candidates, and the true answer
	'''

	from pymongo import MongoClient    # mongodb plugin
	from rdkit import RDLogger
	lg = RDLogger.logger()
	lg.setLevel(4)
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['prediction']
	examples = db['candidates']

	# Define generator
	class Randomizer():
		def __init__(self, seed):
			self.done_ids = []
			np.random.seed(seed)
			if outfile:
				with open(os.path.join(outfile, 'preprocess_candidates_seed.txt'), 'w') as fid:
					fid.write('{}'.format(seed))
		def get_rand(self):
			'''Random WITHOUT replacement'''
			while True:
				doc = examples.find({'found': True, \
					'random': { '$gte': np.random.random()}}).sort('random', 1).limit(1)[0]
				if doc['_id'] in self.done_ids: continue
				self.done_ids.append(doc['_id'])
				yield doc

	if seed == None:
		seed = np.random.randint(10000)
	else:
		seed = int(seed)
	randomizer = Randomizer(seed)
	generator = enumerate(randomizer.get_rand())

	# Initialize (this is not the best way to do this...)
	reaction_strings = []
	reaction_true_onehot = []
	for i, reaction in generator:
		if i == n: break

		strings = [str(reaction['reactant_smiles']) + '>>' + str(x) \
			for x in reaction['product_smiles_candidates']]
		bools = [reaction['product_smiles_true'] in x \
				for x in reaction['product_smiles_candidates']]
		print('rxn. {} : {} true entries'.format(i, sum(bools)))
		reaction_strings.append(
			[x for (y, x) in sorted(zip(bools, strings))]
		)
		reaction_true_onehot.append(
			[y for (y, x) in sorted(zip(bools, strings))]
		)

	return reaction_strings, reaction_true_onehot

def simple_embedding(reaction_strings, reaction_true_onehot, padUpTo = 10000):
	x = []
	y = []

	for i, candidates in enumerate(reaction_strings):
		print('processing reaction {}'.format(i))
		this_rxn = np.vstack(
		 			[rxn_level_descriptors(AllChem.ReactionFromSmarts(z)) for z in candidates]
		 		)
		coo = coo_matrix(pad_sequences(this_rxn.T, maxlen = padUpTo), dtype = np.bool)
		x.append((coo.data, coo.row, coo.col, coo.shape))

	coo = coo_matrix(pad_sequences(reaction_true_onehot, maxlen = padUpTo, dtype = np.bool))
	y = (coo.data, coo.row, coo.col, coo.shape)

	del reaction_strings 
	del reaction_true_onehot

	return x, y

if __name__ == '__main__':
	reaction_strings, reaction_true_onehot = get_candidates(n = 15)
	# for i, example in  enumerate(reaction_true_onehot):
		# print('rxn {}. {} true candidates'.format(i, sum(example)))
		# #print([reaction_strings[i][j] for j, o in enumerate(example) if o])
	x, y = simple_embedding(reaction_strings, reaction_true_onehot, padUpTo = 500)
	# for example in x:
	# 	print(example.shape)

	with open('x_coo.dat', 'wb') as outfile:
		pickle.dump(x, outfile, pickle.HIGHEST_PROTOCOL)
	with open('y_coo.dat', 'wb') as outfile:
		pickle.dump(y, outfile, pickle.HIGHEST_PROTOCOL)
