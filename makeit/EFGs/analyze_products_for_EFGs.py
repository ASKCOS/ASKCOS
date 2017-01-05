import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from pymongo import MongoClient
from bson.objectid import ObjectId
import numpy as np
from EFGs_match import *
from tqdm import tqdm 
import cPickle as pickle
import os

'''This script is meant to test out an EFG-based filtering system
for retrosynthesis. Instead of applying all templates that match, 
we can develop a 'conservative' system by not applying templates in
cases where we have not seen a certain functional group co-present with
the reaction core (in the product).'''

def analyze_products(existence = None, analyzed_rxids = None):
	'''Analyze one transform document'''
	
	# Default values - have not started analysis before
	if existence == None:
		existence = np.array([0.0 for i in range(len(library))])
	if analyzed_rxids == None:
		analyzed_rxids = set()
	total_count = len(analyzed_rxids)
	min_id = max(analyzed_rxids) if analyzed_rxids else 0
	counter = 0

	try:
		for rx_doc in tqdm(REACTION_DB.find({'_id': {'$gt': min_id}}, ['_id', 'RXN_SMILES'],
				no_cursor_timeout = True)):
			if 'RXN_SMILES' not in rx_doc: 
				print('Reaction document does not contain SMILES')
				continue
			if rx_doc['_id'] in analyzed_rxids: continue

			# Create molecule object from the product
			product = Chem.MolFromSmiles(str(rx_doc['RXN_SMILES']).split('>')[-1])
			if not product: 
				print('Could not load product')
				continue 
			product.UpdatePropertyCache()

			# Get EFG signature
			EFG = get_EFGS_matches(product, library = library, exclude = None)
			
			# Update tracker
			existence += np.array(EFG)
			total_count += 1.0
			analyzed_rxids.add(rx_doc['_id'])
			counter += 1

			if counter % 2000 == 0:
				with open('EFG_analyzed_rxids.pickle', 'wb') as fid:
					pickle.dump(analyzed_rxids, fid, protocol = pickle.HIGHEST_PROTOCOL)
				with open('EFG_existence.pickle', 'wb') as fid:
					pickle.dump(existence, fid, protocol = pickle.HIGHEST_PROTOCOL)
				print('Saved!')
	except KeyboardInterrupt:
		print('User requested to stop early')
	except Exception as e:
		print(e)
		print('Something happened? Saving and stopping')

	with open('EFG_analyzed_rxids.pickle', 'wb') as fid:
		pickle.dump(analyzed_rxids, fid, protocol = pickle.HIGHEST_PROTOCOL)
	with open('EFG_existence.pickle', 'wb') as fid:
		pickle.dump(existence, fid, protocol = pickle.HIGHEST_PROTOCOL)

	return existence, total_count

if __name__ == '__main__':


	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['askcos_transforms']
	EFG_DB = db['EFGs']
	library = [doc for doc in EFG_DB.find()]

	# Get retro transforms and instance DB
	db = client['reaxys']
	CHEMICAL_DB = db['chemicals']
	REACTION_DB = db['reactions']


	existence = None
	analyzed_rxids = None
	if os.path.isfile('EFG_analyzed_rxids.pickle') and os.path.isfile('EFG_existence.pickle'):
		with open('EFG_analyzed_rxids.pickle', 'rb') as fid:
			analyzed_rxids = pickle.load(fid)
		with open('EFG_existence.pickle', 'rb') as fid:
			existence = pickle.load(fid)
		print('Reloaded {} previously read RXIDs'.format(len(analyzed_rxids)))

	existence, total_count = analyze_products(existence = existence, analyzed_rxids = analyzed_rxids)
	print(existence)

	# Save to file
	with open('EFG_analysis_products.csv', 'w') as fid:
		fid.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
			'EFG Index', 'Name', 'SMARTS', 'In reduced set?', 'Number of products with group', 
			'Fraction of products with group'
		))

		for i, efg in enumerate(library):
			fid.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
				efg['_id'], efg['name'], efg['SMARTS'], efg['redux'], existence[i], existence[i] / total_count
			))