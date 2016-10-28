from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import sys
import os
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import re
import itertools
from makeit.embedding.reaction_types import get_oxidation_change


if __name__ == '__main__':
	out_folder = os.path.join(os.path.dirname(__file__), 'output')

	# DATABASE
	from pymongo import MongoClient    # mongodb plugin
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaxys']
	TRANSFORM_DB = db['transforms_forward_v1']
	INSTANCE_DB = db['instances']
	CHEMICAL_DB = db['chemicals']
	REACTION_DB = db['reactions']
	REDOX_DB = db['redox']
	done_references = [doc['references'] for doc in REDOX_DB.find({}, ['references'])]
	done_references = list(itertools.chain.from_iterable(done_references))
	max_done_rxid = max([int(rxd_id.split('-')[0]) for rxd_id in done_references])
	print('Largest RXID completed: {}'.format(max_done_rxid))
	made_it_to_new_examples = False

	for reaction in REACTION_DB.find({'_id': {'$gt': max_done_rxid - 1}}, no_cursor_timeout=True):
		# print(reaction)
		# print(str(reaction['RXN_SMILES']))

		reactant_fragments = (str(reaction['RXN_SMILES']).split('>')[0]).split('.')
		product_fragments = (str(reaction['RXN_SMILES']).split('>')[2]).split('.')

		# Largest reactant == one with the most atom mapped atoms
		# largest_reactant = max(reactant_fragments, 
		# 	key = lambda x: len(re.findall('\:[0-9]+\]', x, re.DOTALL)))
		largest_reactant = '.'.join(reactant_fragments) # use all?
		largest_product = max(product_fragments, key = len)

		# print('Primary reactant: {}'.format(largest_reactant))
		# print ('  (with {} mapped atoms'.format(len(re.findall('\:[0-9]+\]', largest_reactant, re.DOTALL))))
		# print('Primary product: {}'.format(largest_product))

		reactant = Chem.MolFromSmiles(largest_reactant)
		product = Chem.MolFromSmiles(largest_product)

		if not reactant:
			print('Could not parse reactant')
			continue
		if not product:
			print('Could not parse product')
			continue

		# Find oxidation change
		delta = get_oxidation_change(reactant, product)

		print('CHANGE: {}'.format(delta))
		if delta == 0:
			print('NEUTRAL! {}'.format(reaction['RXN_SMILES']))
		elif delta < 0:
			print('### REDUCTION ### {}'.format(reaction['RXN_SMILES']))
		elif delta > 0:
			print('### OXIDATION ###  {}'.format(reaction['RXN_SMILES']))

		for i in range(reaction['RX_NVAR']):
			rxd_id = str(reaction['_id']) + '-' + str(i + 1)
			if not made_it_to_new_examples:
				if rxd_id in done_references: 
					continue # save time if already done
				else:
					made_it_to_new_examples = True
			rxd = INSTANCE_DB.find_one({'_id': rxd_id})
			if not rxd: continue

			for xrn in rxd['RXD_RGTXRN']:
				doc = REDOX_DB.find_one({'_id': xrn})
				if not doc:
					new_doc = CHEMICAL_DB.find_one({'_id': xrn}, ['_id', 'SMILES'])
					if not new_doc: continue
					new_doc['references'] = [rxd_id]
					new_doc['count'] = 1.0
					new_doc['ox_power'] = float(delta)
					REDOX_DB.insert(new_doc)
				elif rxd_id not in doc['references']:
					new_power = (delta + doc['ox_power'] * doc['count']) / float(doc['count'] + 1.0)
					REDOX_DB.update_one(
						{'_id': xrn},
						{
							'$set': {
								'ox_power': new_power,
							},
							'$inc': {
								'count': 1,
							},
							'$addToSet': {
								'references': rxd_id,
							},
						}
					)