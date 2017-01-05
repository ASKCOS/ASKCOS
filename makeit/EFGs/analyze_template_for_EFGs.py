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

	reactant_smarts = Chem.MolFromSmarts(str(tform_doc['reaction_smarts'].split('>')[0]))
	coexistence = np.array([0.0 for i in range(len(library))])
	total_count = 0.0

	for rxd_id in tqdm(tform_doc['references']):

		# Create molecule object from the product
		rx_doc = REACTION_DB.find_one({'_id': int(rxd_id.split('-')[0])}, ['RX_PXRN'])
		pxrns = rx_doc['RX_PXRN']
		product_list = []
		for xrn in pxrns:
			chem_doc = CHEMICAL_DB.find_one({'_id': xrn})
			if not chem_doc: continue
			product_list.append(str(chem_doc['SMILES']))
		product = Chem.MolFromSmiles('.'.join(product_list))
		product.UpdatePropertyCache()
		if not product: continue 

		# Figure out which atoms match template
		product_match_ids = product.GetSubstructMatch(reactant_smarts)
		if not product_match_ids:
			print('Could not match reaction back to template?')
			print('Template: {}'.format(str(tform_doc['reaction_smarts'].split('>')[0])))
			print('Reactants: {}'.format(product_list))
			continue 

		# Convert to emol to remove atoms that match templates
		emol = AllChem.RWMol(product)
		for i, idx in enumerate(product_match_ids):
			# Adjust for changed atom IDs from previous atom removals
			offset = sum([idx > old_idx for old_idx in product_match_ids[:i]])
			emol.RemoveAtom(idx - offset)
		mol = emol.GetMol()
		Chem.SanitizeMol(mol,Chem.SanitizeFlags.SANITIZE_PROPERTIES|Chem.SanitizeFlags.SANITIZE_SYMMRINGS)

		# Get EFG signature
		coexistence += np.array(get_EFGS_matches(mol, library = library))
		total_count += 1.0

	coexistence = coexistence / total_count
	return coexistence

if __name__ == '__main__':

	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['askcos_transforms']
	EFG_DB = db['EFGs']
	library = [doc for doc in EFG_DB.find()]

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

	print(coexistence)

	# Save to file
	with open('EFG_analysis_{}.csv'.format(tform_doc['_id']), 'w') as fid:
		fid.write('_id: {}\n'.format(tform_doc['_id']))
		fid.write('SMARTS: {}\n'.format(tform_doc['reaction_smarts']))
		fid.write('number of examples: {}\n'.format(len(tform_doc['references'])))
		fid.write('{}\t{}\t{}\t{}\n'.format(
			'EFG Index', 'Name', 'SMARTS', 'Fraction of examples with co-existing group'
		))

		for i, efg in enumerate(library):
			fid.write('{}\t{}\t{}\t{}\n'.format(
				i, efg['name'], efg['SMARTS'], coexistence[i]
			))