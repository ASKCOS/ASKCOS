from __future__ import print_function
from collections import defaultdict
import rdkit.Chem as Chem 
import rdkit.Chem.AllChem as AllChem
import os
from tqdm import tqdm
import cPickle as pickle
import re

'''

This script is meant to look througuh the Reaxys reaction records and identify
which ones only have partial atom mapping in the products. Unfortunately, even
with RX_SKW: "mapped reaction", some entries are completely unmapped. There is
also the issue of having some reactions with [H] explicit.

Many also have small salt fragments (e.g., counterions) which will be seen as a
contributed-atom, for better or for worse (for worse). 

The script runs locally at first, populating a dictionary with automatic backup, 
before inserting the documents into the database when it's done. For each reaction
with unmapped product fragments, we keep track of which reagents were present in
those RXDs and add them to a list. Ideally, as we aggregate more and more, we can
identify two things:

(1) The most common fragments that are not mapped in the products
- This should be simple additions like [O], [Cl], etc.
- Larger, unpopular fragments will indicate database errors. For example, if
  a 10-atom fragment is unmapped in the producuts, it's safe to say that this
  is a result of Reaxys missing the reactant for that reaction. An example of
  this is RXID 45152

(2) For each fragment, which reagents tended to be present
- This will give us an idea of what reagents are good for, e.g., chlorination
  reactions. 

Note that we do NOT load the structure of the reagents/catalysts/solvents and
determine if the fragment is actually present in those species. 
'''


if __name__ == '__main__':
	out_folder = os.path.join(os.path.dirname(__file__), 'output')
	backup_file = os.path.join(out_folder, 'reagents_contributing_atoms.pickle')

	# DATABASE
	from pymongo import MongoClient    # mongodb plugin
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaxys']
	CHEMICAL_DB = db['chemicals']
	REACTION_DB = db['reactions']
	INSTANCE_DB = db['instances']
	REAGENT_FRAGMENTS_DB = db['reagent_fragments']

	if os.path.isfile(backup_file):
		with open(backup_file, 'rb') as fid:
			as_agent = pickle.load(fid)
		start_at_id = max([max(v['references']) for v in as_agent.values()])
	else:
		as_agent = {} 
		start_at_id = -1

	# Convert to a defaultdict (even if re-initializing)
	as_agent = defaultdict(lambda: {
			'references': set(),
			'reagents': defaultdict(int),
			'catalysts': defaultdict(int),
			'solvents': defaultdict(int),
		}, as_agent)
	print('Starting at ID > {}'.format(start_at_id))
	# first key = the fragment that's contributed
	# reagents = dict mapping xrn to # occurrences

	unmapped = 0
	for i, rx_doc in tqdm(enumerate(
			REACTION_DB.find({ '$and': [
				{'_id': {'$gt': start_at_id}}, {'RX_SKW': 'mapped reaction'}
			]}, ['RX_NVAR', 'RXN_SMILES'], no_cursor_timeout = True))):
		if 'RXN_SMILES' not in rx_doc: continue

		# Save every 50,000
		if i % 50000 == 0:
			with open(backup_file, 'wb') as fid:
				pickle.dump(dict(as_agent), fid, protocol = pickle.HIGHEST_PROTOCOL)
			print('Unmapped reactions: {}/{}'.format(unmapped, i))
		
		prod_smiles = rx_doc['RXN_SMILES'].split('>')[-1]

		# # Make hydrogens implicit first...
		# prod_smiles = prod_smiles.replace('([H])', '')
		# prod_smiles = prod_smiles.replace('[H]', '')

		# # Check for single-atom counterion fragments...
		# prod_smiles = '.'.join([prod_frag for prod_frag in prod_smiles.split('.') 
		# 	if prod_frag not in ['[I-]', '[Cl-]', ['Br-'], ['H+'], ['NH4+']]])

		# Now load into RDKit
		product = Chem.MolFromSmiles(prod_smiles)
		if not product: continue 
		try:
			AllChem.RemoveHs(product)
		except ValueError as e:
			print(e)
			continue

		# Find unmapped atom IDs
		unmapped_ids = [
			a.GetIdx() for a in product.GetAtoms() if not a.HasProp('molAtomMapNumber')
		]
		if not unmapped_ids: continue

		# Define new atom symbols for fragment with atom maps, generalizing fully
		atom_symbols = ['[{}]'.format(a.GetSymbol()) for a in product.GetAtoms()]
		# And bond symbols...
		bond_symbols = ['~' for b in product.GetBonds()]

		extra_reactant_fragment = \
			AllChem.MolFragmentToSmiles(product, unmapped_ids, 
			allHsExplicit = False, isomericSmiles = False, 
			atomSymbols = atom_symbols, bondSymbols = bond_symbols)

		# Treat simultaneous additions as unique
		fragment = '_'.join(sorted(list(set(extra_reactant_fragment.split('.')))))

		# Add RXID reference
		as_agent[fragment]['references'].add(rx_doc['_id'])
					
		# Look up instances of this reaction to see the reagent(s)
		rxd_id_list = ['{}-{}'.format(rx_doc['_id'], j) for j in range(1, int(rx_doc['RX_NVAR']) + 1)]
		for rxd_id in rxd_id_list:
			rxd_doc = INSTANCE_DB.find_one({'_id': rxd_id}, ['RXD_RGTXRN', 'RXD_CATXRN', 'RXD_SOLXRN'])
			if not rxd_doc: continue 

			for xrn in rxd_doc['RXD_RGTXRN']:
				as_agent[fragment]['reagents'][xrn] += 1
			for xrn in rxd_doc['RXD_CATXRN']:
				as_agent[fragment]['catalysts'][xrn] += 1
			for xrn in rxd_doc['RXD_SOLXRN']:
				as_agent[fragment]['solvents'][xrn] += 1

		# print(as_agent)
		# raw_input('Pause...')
		unmapped += 1

	# Turn dictionary into list of docs to insert
	docs = []
	for (key, val) in as_agent.iteritems():
		val['_id'] = key
		val['num_atoms'] = key.count('[')
		val['count'] = len(val['references'])
		doc.append(val)

	REAGENT_FRAGMENTS_DB.insert_many(docs)