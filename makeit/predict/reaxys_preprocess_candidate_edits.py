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
import argparse

FROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data_edits')


def string_or_range_to_float(text):
	try:
		return float(text)
	except Exception as e:
		if '-' in text:
			try:
				return sum([float(x) for x in text.split('-')]) / len(text.split('-'))
			except Exception as e:
				print(e)
		else:
			print(e)
	return None

def get_candidates(candidate_collection, n = 2, seed = None, outfile = '.', shuffle = False, 
	skip = 0, padUpTo = 500, maxEditsPerClass = 5):
	'''
	Pull n example reactions, their candidates, and the true answer
	'''

	from pymongo import MongoClient    # mongodb plugin
	from rdkit import RDLogger
	lg = RDLogger.logger()
	lg.setLevel(4)
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['prediction']
	examples = db[candidate_collection]

	db = client['reaxys']
	INSTANCE_DB = db['instances']
	CHEMICAL_DB = db['chemicals']
	SOLVENT_DB = db['solvents']

	# Define generator
	class Randomizer():
		def __init__(self, seed):
			self.done_ids = []
			self.done_smiles = []
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
					if doc[0]['reactant_smiles'] in self.done_smiles: 
						print('New ID {}, but old reactant SMILES {}'.format(doc[0]['_id'], doc[0]['reactant_smiles']))
						continue
					self.done_ids.append(doc[0]['_id'])
					self.done_smiles.append(doc[0]['reactant_smiles'])
					yield doc[0]
				except KeyboardInterrupt:
					print('Terminated early')
					quit(1)
				except:
					pass

		def get_sequential(self):
			'''Sequential'''
			for doc in examples.find({'found': True}, no_cursor_timeout = True):
				try:
					if not doc: continue 
					if doc['_id'] in self.done_ids: continue
					#if doc['reactant_smiles'] in self.done_smiles: 
					#	print('New ID {}, but old reactant SMILES {}'.format(doc['_id'], doc['reactant_smiles']))
					#	continue
					self.done_ids.append(doc['_id'])
					#self.done_smiles.append(doc['reactant_smiles'])
					yield doc
				except KeyboardInterrupt:
					print('Terminated early')
					quit(1)

	if seed == None:
		seed = np.random.randint(10000)
	else:
		seed = int(seed)
	randomizer = Randomizer(seed)
	if shuffle:
		generator = enumerate(randomizer.get_rand())
	else:
		generator = enumerate(randomizer.get_sequential())

	# Initialize (this is not the best way to do this...)
	reaction_candidate_edits = []
	reaction_candidate_smiles = []
	reaction_true_onehot = []
	reaction_true = []
	reaction_contexts = []

	for i, reaction in generator:
		if i < skip: continue
		if i == skip + n: break

		candidate_smiles = [a for (a, b) in reaction['edit_candidates']]
		candidate_edits =    [b for (a, b) in reaction['edit_candidates']]
		reactant_smiles = reaction['reactant_smiles']
		product_smiles_true = reaction['product_smiles_true']

		# Make sure number of edits is acceptable
		valid_edits = [all([len(e) <= maxEditsPerClass for e in es]) for es in candidate_edits]
		print('Total number of edit candidates: {}'.format(len(valid_edits)))
		print('Total number of valid edit candidates: {}'.format(sum(valid_edits)))
		candidate_smiles = [a for (j, a) in enumerate(candidate_smiles) if valid_edits[j]]
		candidate_edits  = [a for (j, a) in enumerate(candidate_edits) if valid_edits[j]]

		bools = [product_smiles_true == x for x in candidate_smiles]
		print('rxn. {} : {} true entries out of {}'.format(i, sum(bools), len(bools)))
		if sum(bools) > 1:
			print('More than one true? Will take first one')
			# print(reactant_smiles)
			# for (edit, prod) in [(edit, prod) for (boolean, edit, prod) in zip(bools, candidate_edits, candidate_smiles) if boolean]:
			# 	print(prod)
			# 	print(edit)
			# raw_input('Pausing...')
			# continue
			pass
		if sum(bools) == 0:
			print('##### True product not found / filtered out #####')
			continue

		# Sort together and append
		zipsort = sorted(zip(bools, candidate_smiles, candidate_edits))
		zipsort = [[(y, z, x) for (y, z, x) in zipsort if y == 1][0]] + \
				  [(y, z, x) for (y, z, x) in zipsort if y == 0]
		zipsort = zipsort[:padUpTo]

		if sum([y for (y, z, x) in zipsort]) != 1:
			print('New sum true: {}'.format(sum([y for (y, z, x) in zipsort])))
			print('## wrong number of true results?')
			raw_input('Pausing...')

		reaction_candidate_edits.append([x for (y, z, x) in zipsort])
		reaction_true_onehot.append([y for (y, z, x) in zipsort])
		reaction_candidate_smiles.append([z for (y, z, x) in zipsort])
		reaction_true.append(str(reactant_smiles) + '>>' + str(product_smiles_true) + '[{}]'.format(len(zipsort)))

		### Look for conditions
		rxd = INSTANCE_DB.find_one({'_id': reaction['_id']})
		if not rxd:
			print('Could not find RXD with ID {}'.format(reaction['_id']))
			continue
		
		# Temp
		T = string_or_range_to_float(rxd['RXD_T'])
		if not T: T = 20
		# Solvent(s)
		solvent = [0, 0, 0, 0, 0, 0] # c, e, s, a, b, v
		unknown_solvents = []
		for xrn in rxd['RXD_SOLXRN']:
			smiles = str(CHEMICAL_DB.find_one({'_id': xrn})['SMILES'])
			if not smiles: continue 
			mol = Chem.MolFromSmiles(smiles)
			if not mol: continue
			smiles = Chem.MolToSmiles(mol)
			doc = SOLVENT_DB.find_one({'_id': smiles})
			if not doc: 
				unknown_solvents.append(smiles)
				print('Solvent {} not found in DB'.format(smiles))
				continue
			solvent[0] += doc['c']
			solvent[1] += doc['e']
			solvent[2] += doc['s']
			solvent[3] += doc['a']
			solvent[4] += doc['b']
			solvent[5] += doc['v']
		if solvent == [0, 0, 0, 0, 0, 0]:
			doc = SOLVENT_DB.find_one({'_id': 'default'})
			solvent = [doc['c'], doc['e'], doc['s'], doc['a'], doc['b'], doc['v']]
			print('Because all solvents unknown ({}), using default params'.format(', '.join(unknown_solvents)))
		# Reagents/catalysts (as fingerprint, blegh)
		reagent_fp = np.zeros(256) 
		for xrn in rxd['RXD_RGTXRN'] + rxd['RXD_CATXRN']:
			smiles = str(CHEMICAL_DB.find_one({'_id': xrn})['SMILES'])
			if not smiles: continue 
			mol = Chem.MolFromSmiles(smiles)
			if not mol: continue
			reagent_fp += np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits = 256))
		reaction_contexts.append([T] + solvent + list(reagent_fp))

	return reaction_candidate_edits, reaction_true_onehot, reaction_candidate_smiles, reaction_true, reaction_contexts


if __name__ == '__main__':
	n = 5
	padUpTo = 500
	shuffle = False
	skip = 0

	parser = argparse.ArgumentParser()
	parser.add_argument('-n', '--num', type = int, default = 100,
						help = 'Number of candidate sets to read, default 100')
	parser.add_argument('-p', '--padupto', type = int, default = 100,
						help = 'Number of candidates to allow per example, default 100')
	parser.add_argument('-s', '--shuffle', type = int, default = 0,
						help = 'Whether or not to shuffle, default 0')
	parser.add_argument('--skip', type = int, default = 0,
						help = 'How many entries to skip before reading, default 0')
	parser.add_argument('--candidate_collection', type = str, default = 'reaxys_edits_v1',
						help = 'Name of collection within "prediction" db')
	parser.add_argument('--maxeditsperclass', type = int, default = 5, 
						help = 'Maximum number of edits per edit class, default 5')
	args = parser.parse_args()

	n = int(args.num)
	padUpTo = int(args.padupto)
	shuffle = bool(args.shuffle)
	skip = int(args.skip)
	maxEditsPerClass = int(args.maxeditsperclass)

	reaction_candidate_edits, reaction_true_onehot, reaction_candidate_smiles, reaction_true, reaction_contexts = \
			get_candidates(args.candidate_collection, n = n, shuffle = shuffle, skip = skip, 
				padUpTo = padUpTo, maxEditsPerClass = maxEditsPerClass)

	rounded_time = '{}-{}'.format(skip, skip + n - 1)

	with open(os.path.join(FROOT, '{}_candidate_edits.pickle'.format(rounded_time)), 'wb') as outfile:
		pickle.dump(reaction_candidate_edits, outfile, pickle.HIGHEST_PROTOCOL)
	with open(os.path.join(FROOT, '{}_candidate_bools.pickle'.format(rounded_time)), 'wb') as outfile:
		pickle.dump(reaction_true_onehot, outfile, pickle.HIGHEST_PROTOCOL)
	with open(os.path.join(FROOT, '{}_candidate_smiles.pickle'.format(rounded_time)), 'wb') as outfile:
		pickle.dump(reaction_candidate_smiles, outfile, pickle.HIGHEST_PROTOCOL)
	with open(os.path.join(FROOT, '{}_reaction_string.pickle'.format(rounded_time)), 'wb') as outfile:
		pickle.dump(reaction_true, outfile, pickle.HIGHEST_PROTOCOL)
	with open(os.path.join(FROOT, '{}_contexts.pickle'.format(rounded_time)), 'wb') as outfile:
		pickle.dump(reaction_contexts, outfile, pickle.HIGHEST_PROTOCOL)
