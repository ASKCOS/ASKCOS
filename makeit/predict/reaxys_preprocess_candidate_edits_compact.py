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
from makeit.embedding.descriptors import edits_to_vectors

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


def get_candidates(candidate_collection, outfile = '.', n_max = 500):
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
		def __init__(self):
			self.done_ids = []
			self.done_smiles = []

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

	randomizer = Randomizer()
	generator = enumerate(randomizer.get_sequential())

	i = 0
	for j, reaction in generator:
		
		candidate_smiles = [a for (a, b) in reaction['edit_candidates']]
		candidate_edits =    [b for (a, b) in reaction['edit_candidates']]
		reactant_smiles = reaction['reactant_smiles']
		product_smiles_true = reaction['product_smiles_true']

		reactants_check = Chem.MolFromSmiles(str(reactant_smiles))
		if not reactants_check:
			continue

		bools = [product_smiles_true == x for x in candidate_smiles]
		# print('rxn. {} : {} true entries out of {}'.format(i, sum(bools), len(bools)))
		if sum(bools) > 1:
			# print('More than one true? Will take first one')
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

		#print('Total number of edit candidates: {} ({} valid)'.format(len(valid_edits), sum(valid_edits)))
		
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
			chem_doc = CHEMICAL_DB.find_one({'_id': xrn}, ['SMILES'])
			if not chem_doc: continue
			smiles = str(chem_doc['SMILES'])
			if not smiles: continue 
			mol = Chem.MolFromSmiles(smiles)
			if not mol: continue
			smiles = Chem.MolToSmiles(mol)
			doc = SOLVENT_DB.find_one({'_id': smiles})
			context_info += smiles + ','
			if not doc: 
				unknown_solvents.append(smiles)
				#print('Solvent {} not found in DB'.format(smiles))
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

		reaction_candidate_edits_compact = [
				';'.join([
					','.join(x[0]),
					','.join(x[1]),
					','.join(['%s-%s-%s' % tuple(blost) for blost in x[2]]),
					','.join(['%s-%s-%s' % tuple(bgain) for bgain in x[3]]),
				])
			for (y, z, x) in zipsort]

		# legend = {
		# 	'candidate_edits_compact': 0,
		# 	'atom_desc_dict': 1,
		# 	'T': 2,
		# 	'solvent': 3,
		# 	'reagent': 4,
		# 	'yield': 5,
		# 	'rxdid': 6,
		# 	'reaction_true_onehot': 7,
		# 	'candidate_smiles': 8,
		# 	'reaction_true': 9,
		# }

		if not countonly:
			pickle.dump((
				reaction_candidate_edits_compact,
				edits_to_vectors([], reactants_check, return_atom_desc_dict = True),
				T,
				solvent,
				reagent_fp,
				rxd['RXD_NYD'],
				[y for (y, z, x) in zipsort],
			), fid_data, pickle.HIGHEST_PROTOCOL)

			pickle.dump((
				str(rxd['_id']),
				str(reactant_smiles) + '>' + str(context_info) + '>' + str(product_smiles_true) + '[{}]'.format(len(zipsort)),
				[z for (y, z, x) in zipsort],
				reaction_candidate_edits_compact,
			), fid_labels, pickle.HIGHEST_PROTOCOL)

		i += 1

		if i % 100 == 0:
			print('Completed {}/{}'.format(i, n_max))

		if i == n_max:
			print('Finished the requested {} examples'.format(n_max))
			break


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-n', '--num', type = int, default = 75,
					help = 'Maximum number of examples to save')
	parser.add_argument('--candidate_collection', type = str, default = 'reaxys_edits_v1',
						help = 'Name of collection within "prediction" db')
	parser.add_argument('--complete_only', type = str, default = 'y', 
						help = 'Only use complete examples, including recognized solvent, default y')
	parser.add_argument('--tag', type = str, default = 'reaxys',
						help = 'Data tag, default = reaxys')
	parser.add_argument('--count', type = int, default = 0, help = 'Only count, dont actually save?')
	args = parser.parse_args()

	n = int(args.num)
	complete_only = args.complete_only in ['y', 'Y', 'T', 't', 'true', '1']
	tag = str(args.tag)
	countonly = bool(args.count)
	print('Only using complete records')

	legend_data = {
		'candidate_edits_compact': 0,
		'atom_desc_dict': 1,
		'T': 2,
		'solvent': 3,
		'reagent': 4,
		'yield': 5,
		'reaction_true_onehot': 6,
		'N_examples': n,
	}

	legend_labels = {
		'rxdid': 0,
		'reaction_true': 1,
		'candidate_smiles': 2,
		'candidate_edits_compact': 3,
		'N_examples': n,
	}

	# Open pickle file
	if not countonly:
		fid_data = open(os.path.join(FROOT, '{}_data.pickle'.format(tag)), 'wb')
		fid_labels = open(os.path.join(FROOT, '{}_labels.pickle'.format(tag)), 'wb')

		# First entry is the legend
		pickle.dump(legend_data, fid_data, pickle.HIGHEST_PROTOCOL)
		pickle.dump(legend_labels, fid_labels, pickle.HIGHEST_PROTOCOL)

	# Now get candidates
	get_candidates(args.candidate_collection, n_max = n)

	if not countonly:
		fid_data.close()
		fid_labels.close()
