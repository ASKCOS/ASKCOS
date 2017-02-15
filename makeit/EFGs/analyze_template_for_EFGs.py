import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from pymongo import MongoClient
from bson.objectid import ObjectId
import numpy as np
from EFGs_match import *
from tqdm import tqdm 

'''This script is meant to test out an EFG-based filtering system
for retrosynthesis. Instead of applying all templates that match, 
we can develop a 'conservative' system by not applying templates in
cases where we have not seen a certain functional group co-present with
the reaction core (in the product).'''

def analyze_template(tform_doc):
	'''Analyze one transform document'''

	reactant_smarts = str(tform_doc['reaction_smarts'].split('>')[0])
	coexistence = np.array([0.0 for i in range(len(library))])
	total_count = 0.0

	# TEMP YIELD TEST
	coexistence_yields = [[] for i in range(len(library))]
	total_yield_hasYield = 0.0
	total_count_hasYield = 0.0

	for rxd_id in tqdm(tform_doc['references']):

		# Create molecule object from the product
		rx_doc = REACTION_DB.find_one({'_id': int(rxd_id.split('-')[0])})#, ['RX_PXRN', 'RXN_SMILES'])
		pxrns = rx_doc['RX_PXRN']
		product_list = []
		for xrn in pxrns:
			chem_doc = CHEMICAL_DB.find_one({'_id': xrn})
			if not chem_doc: continue
			product_list.append(str(chem_doc['SMILES']))
		product = Chem.MolFromSmiles('.'.join(product_list))
		product.UpdatePropertyCache()
		if not product: continue 

		# Get EFG signature
		EFG = get_EFGS_matches(product, library = library, exclude = reactant_smarts)
		coexistence += np.array(EFG)
		total_count += 1.0

		# # TEMPORARY FOR DEBUGGING
		# TEST_INDEX = 210
		# if EFG[TEST_INDEX] == 1:
		# 	print('Has a {}'.format(library[TEST_INDEX]['name']))
		# 	print(product_list)
		# 	print(rx_doc)
		# 	raw_input('Pause..')

		# TEMP YIELD TEST
		rxd_doc = INSTANCE_DB.find_one({'_id': rxd_id})
		if 'RXD_NYD' not in rxd_doc: continue 
		for i, exists in enumerate(EFG):
			if exists: coexistence_yields[i].append(rxd_doc['RXD_NYD'])
		total_yield_hasYield += rxd_doc['RXD_NYD']
		total_count_hasYield += 1.0

	np.seterr(divide='ignore', invalid='ignore')
	coexistence = coexistence / total_count
	coexistence_yields_avg = [np.mean(yields) for yields in coexistence_yields]
	coexistence_yields_std = [np.std(yields) for yields in coexistence_yields]

	# Save to file
	with open('{}_analysis_{}.csv'.format(label, tform_doc['_id']), 'w') as fid:
		fid.write('_id: {}\n'.format(tform_doc['_id']))
		fid.write('SMARTS: {}\n'.format(tform_doc['reaction_smarts']))
		fid.write('number of examples: {}\n'.format(len(tform_doc['references'])))
		fid.write('{}\t{}\t{}\t{}\t{}\n'.format(
			'EFG Index', 'Name', 'SMARTS', 'In reduced set?', 'Fraction of examples with co-existing group'
		))

		for i, efg in enumerate(library):
			fid.write('{}\t{}\t{}\t{}\t{}\n'.format(
				efg['_id'], efg['name'], efg['SMARTS'], efg['redux'], coexistence[i]
			))

	# TEMP YIELDS ANALYSIS -  Save to file
	with open('{}_analysis_{}.csv'.format(label, tform_doc['_id']), 'w') as fid:
		fid.write('_id: {}\n'.format(tform_doc['_id']))
		fid.write('SMARTS: {}\n'.format(tform_doc['reaction_smarts']))
		fid.write('number of examples: {}\n'.format(len(tform_doc['references'])))
		fid.write('average yield: {}\n'.format(total_yield_hasYield / total_count_hasYield))
		fid.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
			'EFG Index', 'Name', 'SMARTS', 'In reduced set?', 'Fraction of examples with co-existing group', 
			'Avg yield when group present', 'Stdev of yield when present', 'Num ex with yield'
		))

		for i, efg in enumerate(library):
			fid.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
				efg['_id'], efg['name'], efg['SMARTS'], efg['redux'], coexistence[i], 
				coexistence_yields_avg[i], coexistence_yields_std[i], len(coexistence_yields[i])
			))

	return coexistence

if __name__ == '__main__':

	USE_REDUCED_EFGs = False

	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['askcos_transforms']
	EFG_DB = db['EFGs']
	if USE_REDUCED_EFGs:
		db_filter = {'redux': True}
		label = 'EFG_redux'
	else:
		db_filter = {}
		label = 'EFG'
	library = [doc for doc in EFG_DB.find(db_filter)]

	# Get retro transforms and instance DB
	db = client['reaxys']
	TRANSFORM_DB = db['transforms_retro_v2']
	INSTANCE_DB = db['instances']
	CHEMICAL_DB = db['chemicals']
	REACTION_DB = db['reactions']

	# Get ONE example for now
	tform_doc = TRANSFORM_DB.find_one({'count': {'$gt': 1000, '$lt': 2000}})
	tform_doc = TRANSFORM_DB.find_one({'_id': ObjectId('57c829b5fbff5064589088d0')})
	if not tform_doc: quit(1)
	print('Template SMARTS: {}'.format(str(tform_doc['reaction_smarts'])))
	print('id: {}'.format(tform_doc['_id']))
	coexistence = analyze_template(tform_doc)
