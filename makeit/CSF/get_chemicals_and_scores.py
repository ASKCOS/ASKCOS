# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import numpy as np
import os
import sys
import cPickle as pickle
from collections import defaultdict
from tqdm import tqdm


BUYABLE_BASE_SCORE = 0
REACTION_COST = 100
UPPER_LIMIT_SCORE = np.inf

REACTION_LIMIT = 5000000
MAX_CHEMICALS = 100

def func_UPPER_LIMIT_SCORE(x = 'optional'):
	return UPPER_LIMIT_SCORE

if __name__ == '__main__':

	FROOT = os.path.dirname(os.path.realpath(__file__))

	from pymongo import MongoClient    # mongodb plugin
	from rdkit import RDLogger
	lg = RDLogger.logger()
	lg.setLevel(4)
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaxys']
	CHEMICAL_DB = db['chemicals']
	REACTION_DB = db['reactions']

	xrn_to_value = defaultdict(func_UPPER_LIMIT_SCORE)
	xrn_to_smiles = {}
	xrn_to_fp = {}

	if os.path.isfile(os.path.join(FROOT, 'initial_chemicals.pickle')):
		print('Reloading previously saved chemicals')
		with open(os.path.join(FROOT, 'initial_chemicals.pickle'), 'rb') as fid:
			xrn_to_value  = pickle.load(fid)
	else:
		# Populate dictionary of chemicals with known buyables
		print('Getting chemicals')
		for i, chem_doc in enumerate(CHEMICAL_DB.find({'buyable_id': {'$gt': -2}}, ['_id'])):

			if i % 10 == 0: print(i)
			xrn_to_value[chem_doc['_id']] = BUYABLE_BASE_SCORE
				
		print('Found {} buyable chemicals from DB'.format(len(xrn_to_value)))

		with open(os.path.join(FROOT, 'initial_chemicals.pickle'), 'wb') as fid:
			pickle.dump(xrn_to_value, fid, pickle.HIGHEST_PROTOCOL)

	if os.path.isfile(os.path.join(FROOT, 'firstpass_chemicals.pickle')):
		print('Reloading previously saved chemicals (after going through reactions)')
		with open(os.path.join(FROOT, 'firstpass_chemicals.pickle'), 'rb') as fid:
			xrn_to_value  = pickle.load(fid)

	else:

		def go_through_once():
			# Look through the reactions
			for rxn_doc in tqdm(REACTION_DB.find({}, ['RX_RXRN', 'RX_PXRN']).limit(REACTION_LIMIT)):
				if len(rxn_doc['RX_RXRN']) == 0: continue 
				if len(rxn_doc['RX_PXRN']) == 0: continue 
				# What is the most expensive reactant required for this reaction?
				max_reactant_score = BUYABLE_BASE_SCORE
				for rx_rxrn in rxn_doc['RX_RXRN']:
					max_reactant_score = max(xrn_to_value[rx_rxrn], max_reactant_score)
				# If one of the reactants isn't buyable, nothing to do
				if max_reactant_score == UPPER_LIMIT_SCORE: continue 
				for rx_pxrn in rxn_doc['RX_PXRN']:
					xrn_to_value[rx_pxrn] = min(xrn_to_value[rx_pxrn], max_reactant_score + REACTION_COST)


						# Evaluate scores
			score_histogram = defaultdict(int)
			for (xrn, val) in xrn_to_value.iteritems():
				score_histogram[val] += 1
			print(score_histogram)

		print('Getting reactions (first time)')
		go_through_once()
		print('Getting reactions (second time)')
		go_through_once()
		print('Getting reactions (third time)')
		go_through_once()
	
		print('Saving reactions (no fingerprints)')
		with open(os.path.join(FROOT, 'firstpass_chemicals.pickle'), 'wb') as fid:
			pickle.dump(xrn_to_value, fid, pickle.HIGHEST_PROTOCOL)

	# Now go back and get fingerprints for other chemicals
	print('Getting fingerprints for chemicals and saving')
	from random import shuffle
	xrns = xrn_to_value.keys()
	shuffle(xrns)
	with open(os.path.join(FROOT, 'data.pickle'), 'wb') as outfid:
		counter = 0
		for i, xrn in tqdm(enumerate(xrns)):
			if xrn_to_value[xrn] == np.inf: continue

			if counter == MAX_CHEMICALS: break

			chem_doc = CHEMICAL_DB.find_one({'_id': xrn}, ['SMILES', 'M2_fp1024'])
			if not chem_doc: continue
			if 'M2_fp1024' not in chem_doc: continue
			if 'SMILES' not in chem_doc: continue

			smiles = chem_doc['SMILES']
			fp = np.unpackbits(pickle.loads(chem_doc['M2_fp1024']))

			pickle.dump((xrn, smiles, fp, xrn_to_value[xrn]), outfid)
			counter += 1
