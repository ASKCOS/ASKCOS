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

FROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data_edits_reaxys')

def string_or_range_to_float(text):
	try:
		return float(text)
	except Exception as e:
		x = [z for z in text.strip().split('-') if z not in [u'', u' ']]
		if text.count('-') == 1: # 20 - 30
			try:
				return (float(x[0]) + float(x[1])) / 2.0
			except Exception as e:
				print('Could not convert {}, {}'.format(text, x))
				#print(e)
		elif text.count('-') == 2: # -20 - 0
			try:
				return (-float(x[0]) + float(x[1])) / 2.0
			except Exception as e:
				print('Could not convert {}, {}'.format(text, x))
				#print(e)
		elif text.count('-') == 3: # -20 - -10
			try:
				return (-float(x[0]) - float(x[1])) / 2.0
			except Exception as e:
				print('Could not convert {}, {}'.format(text, x))
				#print(e)
		else:
			print('Could not convert {}'.format(text))
			print(e)
	return np.nan


def get_candidates(candidate_collection, seed = 0, outfile = '.', shuffle = False, 
	skip = 0, padUpTo = 500, maxEditsPerClass = 5, n_max = 500):
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
			for doc in examples.find({'found': True}, no_cursor_timeout = True).sort('_id', 1):
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

	randomizer = Randomizer(seed)
	if shuffle:
		generator = enumerate(randomizer.get_rand())
	else:
		generator = enumerate(randomizer.get_sequential())

	### Initialize 
	from makeit.embedding.descriptors import edits_to_vectors
	mol = Chem.MolFromSmiles('[C:1][C:2]')
	(a, _, b, _) = edits_to_vectors((['1'],[],[('1','2',1.0)],[]), mol)
	F_atom = len(a[0])
	F_bond = len(b[0])


	import h5py
	h5f = h5py.File(os.path.join(FROOT, 'data_edits_reaxys.h5'), 'w')

	asciistring = h5py.special_dtype(vlen = bytes) # text
	regularfloat = np.float32

	h5f.create_dataset('x_h_lost', (n_max, padUpTo, maxEditsPerClass, F_atom), dtype=regularfloat)
	h5f.create_dataset('x_h_gain', (n_max, padUpTo, maxEditsPerClass, F_atom), dtype=regularfloat)
	h5f.create_dataset('x_bond_lost', (n_max, padUpTo, maxEditsPerClass, F_bond), dtype=regularfloat)
	h5f.create_dataset('x_bond_gain', (n_max, padUpTo, maxEditsPerClass, F_bond), dtype=regularfloat)
	h5f.create_dataset('T', (n_max,), dtype=regularfloat)
	h5f.create_dataset('solvent', (n_max, 6), dtype=regularfloat)
	h5f.create_dataset('reagent', (n_max, 256), dtype=regularfloat)
	h5f.create_dataset('yield', (n_max,), dtype=regularfloat)

	h5f.create_dataset('reaction_true_onehot', (n_max,padUpTo), dtype=bool)
	h5f.create_dataset('rxdid', (n_max,), dtype=asciistring)
	h5f.create_dataset('candidate_smiles', (n_max,padUpTo), dtype=asciistring)
	h5f.create_dataset('candidate_edits', (n_max,padUpTo), dtype=asciistring)
	h5f.create_dataset('reaction_true', (n_max,), dtype=asciistring)

	i = 0
	for j, reaction in generator:
		
		candidate_smiles = [a for (a, b) in reaction['edit_candidates']]
		candidate_edits =    [b for (a, b) in reaction['edit_candidates']]
		reactant_smiles = reaction['reactant_smiles']
		product_smiles_true = reaction['product_smiles_true']

		reactants_check = Chem.MolFromSmiles(str(reactant_smiles))
		if not reactants_check:
			continue

		# Make sure number of edits is acceptable
		valid_edits = [all([len(e) <= maxEditsPerClass for e in es]) for es in candidate_edits]
		candidate_smiles = [a for (j, a) in enumerate(candidate_smiles) if valid_edits[j]]
		candidate_edits  = [a for (j, a) in enumerate(candidate_edits) if valid_edits[j]]

		bools = [product_smiles_true == x for x in candidate_smiles]
		# print('rxn. {} : {} true entries out of {}'.format(i, sum(bools), len(bools)))
		if sum(bools) > 1:
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

		### Look for conditions
		context_info = ''
		rxd = INSTANCE_DB.find_one({'_id': reaction['_id']})
		if not rxd:
			print('Could not find RXD with ID {}'.format(reaction['_id']))
			raise ValueError('Candidate reaction source not found?')
		if complete_only and 'complete' not in rxd:
			continue

		print('Total number of edit candidates: {} ({} valid)'.format(len(valid_edits), sum(valid_edits)))
		
		# Temp
		T = string_or_range_to_float(rxd['RXD_T'])
		if not T: 
			T = 20
			if complete_only: continue # skip if T was unparseable
		if T == -1: 
			T = 20
			if complete_only: continue # skip if T not recorded

		# Solvent(s)
		solvent = [0, 0, 0, 0, 0, 0] # c, e, s, a, b, v
		unknown_solvents = []
		context_info += 'solv:'
		for xrn in rxd['RXD_SOLXRN']:
			try:
				smiles = str(CHEMICAL_DB.find_one({'_id': xrn})['SMILES'])
			except:
				continue
			if not smiles: continue 
			mol = Chem.MolFromSmiles(smiles)
			if not mol: continue
			smiles = Chem.MolToSmiles(mol)
			doc = SOLVENT_DB.find_one({'_id': smiles})
			context_info += smiles + ','
			if not doc: 
				unknown_solvents.append(smiles)
				print('Solvent {} not found in DB'.format(smiles))
				context_info = context_info[:-1] + '?,' # add question mark to denote unfound
				continue
			solvent[0] += doc['c']
			solvent[1] += doc['e']
			solvent[2] += doc['s']
			solvent[3] += doc['a']
			solvent[4] += doc['b']
			solvent[5] += doc['v']
		if solvent == [0, 0, 0, 0, 0, 0]:
			if complete_only: continue # if solvent not parameterized, skip
			doc = SOLVENT_DB.find_one({'_id': 'default'})
			solvent = [doc['c'], doc['e'], doc['s'], doc['a'], doc['b'], doc['v']]
			print('Because all solvents unknown ({}), using default params'.format(', '.join(unknown_solvents)))
		# Reagents/catalysts (as fingerprint, blegh)
		reagent_fp = np.zeros(256) 
		context_info += 'rgt:'
		for xrn in rxd['RXD_RGTXRN'] + rxd['RXD_CATXRN']:
			doc = CHEMICAL_DB.find_one({'_id': xrn})
			if not doc:
				print('########## COULD NOT FIND REAGENT {} ###########'.format(xrn))
				continue
			smiles = str(doc['SMILES'])
			if not smiles: continue 
			mol = Chem.MolFromSmiles(smiles)
			if not mol: continue
			context_info += smiles + ','
			reagent_fp += np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits = 256))
		context_info += 'T:{}'.format(T)
		if complete_only: # should have info about time and yield
			context_info += ',t:' + str(rxd['RXD_TIM']) + 'min,y:' + str(rxd['RXD_NYD']) + '%'
		#print(context_info)


		### NOW CONVERT THE DATA TO NUMERICAL VALUES ###
		x_h_lost = np.zeros((padUpTo, maxEditsPerClass, F_atom))
		x_h_gain = np.zeros((padUpTo, maxEditsPerClass, F_atom))
		x_bond_lost = np.zeros((padUpTo, maxEditsPerClass, F_bond))
		x_bond_gain = np.zeros((padUpTo, maxEditsPerClass, F_bond))

		for (c, (z, y, edits)) in enumerate(zipsort):
			if any([len(edit) > 5 for edit in edits]):
				#print('Edit counts: {}'.format([len(edit) for edit in edits]))
				#print('skipping')
				continue
			try:
				edit_h_lost_vec, edit_h_gain_vec, \
					edit_bond_lost_vec, edit_bond_gain_vec = \
					edits_to_vectors(edits, reactants_check)
			except KeyError as e:
				print(e)
				print('Skipping this edit')
				continue
			for (e, edit_h_lost) in enumerate(edit_h_lost_vec):
				x_h_lost[c, e, :] = edit_h_lost
			for (e, edit_h_gain) in enumerate(edit_h_gain_vec):
				x_h_gain[c, e, :] = edit_h_gain
			for (e, edit_bond_lost) in enumerate(edit_bond_lost_vec):
				x_bond_lost[c, e, :] = edit_bond_lost
			for (e, edit_bond_gain) in enumerate(edit_bond_gain_vec):
				x_bond_gain[c, e, :] = edit_bond_gain

		# Get rid of NaNs
		x_h_lost[np.isnan(x_h_lost)] = 0.0
		x_h_gain[np.isnan(x_h_gain)] = 0.0
		x_bond_lost[np.isnan(x_bond_lost)] = 0.0
		x_bond_gain[np.isnan(x_bond_gain)] = 0.0
		x_h_lost[np.isinf(x_h_lost)] = 0.0
		x_h_gain[np.isinf(x_h_gain)] = 0.0
		x_bond_lost[np.isinf(x_bond_lost)] = 0.0
		x_bond_gain[np.isinf(x_bond_gain)] = 0.0

		h5f['x_h_lost'][i] = x_h_lost
		h5f['x_h_gain'][i] = x_h_gain
		h5f['x_bond_lost'][i] = x_bond_lost
		h5f['x_bond_gain'][i] = x_bond_gain
		h5f['T'][i] = T
		h5f['solvent'][i] = np.array(solvent)
		h5f['reagent'][i] = reagent_fp
		h5f['yield'][i] = rxd['RXD_NYD']
		h5f['rxdid'][i] = rxd['_id']

		for c, (y, z, x) in enumerate(zipsort):
			h5f['reaction_true_onehot'][i,c] = bool(y)
			h5f['candidate_smiles'][i,c] = bytes(z)
			h5f['candidate_edits'][i,c] = bytes(x)
		for c_remaining in range(c + 1, padUpTo):
			print(c_remaining)
			h5f['reaction_true_onehot'][i,c_remaining] = False
			h5f['candidate_smiles'][i,c_remaining] = bytes('n.r.')
			h5f['candidate_edits'][i,c_remaining] = bytes('n.r.')
		reaction_true =  str(reactant_smiles) + '>' + str(context_info) + \
			';Nc:{}>'.format(len(zipsort)) + str(product_smiles_true)
		h5f['reaction_true'][i] = bytes(reaction_true)

		i += 1
		if i % 10 == 0:
			print('Logged {} total examples so far'.format(i))

		if i == n_max:
			print('Finished the requested {} examples'.format(n_max))
			h5f.close()
			break


if __name__ == '__main__':
	padUpTo = 500
	shuffle = False
	skip = 0

	parser = argparse.ArgumentParser()
	parser.add_argument('--max', type = int, default = 10000,
					help = 'Maximum number of examples to save')
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
	parser.add_argument('--complete_only', type = str, default = 'y', 
						help = 'Only use complete examples, including recognized solvent, default y')
	args = parser.parse_args()

	padUpTo = int(args.padupto)
	shuffle = bool(args.shuffle)
	skip = int(args.skip)
	maxEditsPerClass = int(args.maxeditsperclass)
	complete_only = args.complete_only in ['y', 'Y', 'T', 't', 'true', '1']
	print('Only using complete records')

	get_candidates(args.candidate_collection, shuffle = shuffle, skip = skip, seed = 0,
				padUpTo = padUpTo, maxEditsPerClass = maxEditsPerClass, n_max = int(args.max))
