# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import argparse
import numpy as np     	      	   # for simple calculations
import os                          # for saving
import sys
import rdkit.Chem as Chem
import makeit.retro.transformer as transformer 
from pymongo import MongoClient    # mongodb plugin
from makeit.predict.summarize_reaction_outcome import summarize_reaction_outcome
from makeit.embedding.descriptors import edits_to_vectors, oneHotVector # for testing
import re
import time
from tqdm import tqdm
from multiprocessing import Process, Lock, Queue
from Queue import Empty as QueueEmpty

'''
This script generates candidates and logs their "candidate_edits" in a prediction 
database. Multiprocessing is used to help speed things up.

Unfortunately, a memory leak in RDKit means that this script will need to be interrupted
and restarted many times...

v2 does not allow solvents or catalysts to react, which should give huge time savings
'''


def process_one(rx, rxd, lock_pymongo):

	# LOGGING
	start_time = time.time()

	# print('#########')
	# print('## RXN {}'.format(rxd['_id']))
	# print('#########')
	print('[RXD ID {}] Starting...'.format(rxd['_id']))

	rxn_smiles = rx['RXN_SMILES']
	reactants = rxn_smiles.split('>')[0]
	products = rxn_smiles.split('>')[2]

	# Add ONLY reagents into the reactant pool
	extra_smiles = ''
	XRNs = rxd['RXD_RGTXRN'] #+ rxd['RXD_SOLXRN'] + rxd['RXD_CATXRN']
	if XRNs:
		docs = [CHEMICAL_DB.find_one({'_id': xrn}) for xrn in XRNs]
		extra_smiles = [doc['SMILES'] for doc in docs if doc]
		extra_smiles = '.'.join([extra_smile for extra_smile in extra_smiles if extra_smile])
	if extra_smiles: reactants += '.' + extra_smiles

	#print('"REACTANTS(+REAGENTS+ETC)": {}'.format(reactants))
	reactants = Chem.MolFromSmiles(reactants)
	products = Chem.MolFromSmiles(products)

	if (not reactants) or (not products):
		print('Could not load reactants/products?')
		print('ID: {}'.format(rxd['_id']))
		return

	[a.ClearProp('molAtomMapNumber') for a in reactants.GetAtoms()] # remove atom mapping
	[a.ClearProp('molAtomMapNumber') for a in products.GetAtoms()] # remove atom mapping

	major_product_smiles = max(Chem.MolToSmiles(products, isomericSmiles = USE_STEREOCHEMISTRY).split('.'), key = len)
	if major_product_smiles in Chem.MolToSmiles(reactants, isomericSmiles = USE_STEREOCHEMISTRY).split('.'):
		#print('WARNING: recorded product already in reactants?')
		return

	#print('REACTANTS: {}'.format(Chem.MolToSmiles(reactants, isomericSmiles = USE_STEREOCHEMISTRY)))
	#print('PRODUCTS:  {}'.format(Chem.MolToSmiles(products, isomericSmiles = USE_STEREOCHEMISTRY)))

	n_reactant_atoms = len(reactants.GetAtoms())
	#print('Number of reactant atoms: {}'.format(n_reactant_atoms))
	if n_reactant_atoms > 80:
		#print('Skipping huge molecule! N_reactant_atoms = {}'.format(n_reactant_atoms))
		return

	# Add new atom map numbers
	[a.SetProp('molAtomMapNumber', str(i+1)) for (i, a) in enumerate(reactants.GetAtoms())]
	
	found_true = False
	candidate_edits = [] # list of tuples of (product smiles, edits required)
	for template in Transformer.templates:
		# Perform transformation
		try:
			outcomes = template['rxn_f'].RunReactants([reactants])
		except Exception as e:
			if v: print(e)
			continue
		if not outcomes: continue # no match
		for j, outcome in enumerate(outcomes):
			outcome = outcome[0] # all products represented as single mol by transforms
			try:
				outcome.UpdatePropertyCache()
				Chem.SanitizeMol(outcome)
				# Figure out what atom mapping was
				# note: this requires a special version of RDKit!!
				[a.SetProp('molAtomMapNumber', a.GetProp('old_molAtomMapNumber')) \
					for (i, a) in enumerate(outcome.GetAtoms()) \
					if 'old_molAtomMapNumber' in a.GetPropsAsDict()]
			except Exception as e:
				#print(e)
				continue
			if v: print('Outcome SMILES: {}'.format(Chem.MolToSmiles(outcome)))

			# Reduce to largest (longest) product only?
			candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY)

			candidate_smiles = max(candidate_smiles.split('.'), key = len)
			outcome = Chem.MolFromSmiles(candidate_smiles)
			
			# Find what edits were made
			edits = summarize_reaction_outcome(reactants, outcome)
			if v: print(edits)

			# Remove mapping before matching
			[x.ClearProp('molAtomMapNumber') for x in outcome.GetAtoms() if x.HasProp('molAtomMapNumber')] # remove atom mapping from outcome
			candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY)
			if candidate_smiles == major_product_smiles:
				if v: print('Matched true [{}]'.format(major_product_smiles))
				found_true = True

			# Add to ongoing list
			if (candidate_smiles, edits) not in candidate_edits:
				candidate_edits.append((candidate_smiles, edits))

	# Prepare doc and insert
	if found_true: print('[RXD ID {}] Found true product'.format(rxd['_id']))
	doc = {
		'_id': rxd['_id'],
		'reactant_smiles': Chem.MolToSmiles(reactants, isomericSmiles = USE_STEREOCHEMISTRY),
		'edit_candidates': candidate_edits,
		'product_smiles_true': major_product_smiles,
		'found': found_true,
		'num_candidates': len(candidate_edits),
	}

	lock_pymongo.acquire()
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['prediction']
	candidates = db[candidate_collection]
	res = candidates.insert(doc)
	lock_pymongo.release()
	print('[RXD ID {}] Inserted candidates into db'.format(doc['_id']))

	# LOGGING
	end_time = time.time()
	print('[RXD ID {}] Time: {}'.format(rxd['_id'], end_time - start_time))
	print('[RXD ID {}] Unique edit sets using longest prod: {}'.format(rxd['_id'], len(candidate_edits)))

	return

def process_forever(queue, lock_pymongo):
	while True:
		try:
			(rx, rxd) = queue.get()
			process_one(rx, rxd, lock_pymongo)
		except QueueEmpty:
			time.sleep(1)

def rxd_generator():
	'''Return (rx, rxd) tuples that have not been processed'''

	print('Looking for done_ids for the first time...')
	done_ids = set(
		[doc['_id'] for doc in CANDIDATE_DB.find({}, ['_id'], no_cursor_timeout = True)]
	)
	print('...done')

	for rx in REACTION_DB.find({'RX_SKW': 'mapped reaction'}, ['_id', 'RXN_SMILES', 'RX_NVAR'], 
			no_cursor_timeout = True).sort('_id', 1):
		if not rx: break 
		for i in range(1, rx['RX_NVAR'] + 1):
			rxd_id = '{}-{}'.format(rx['_id'], i)
			if rxd_id in done_ids: continue
			rxd = INSTANCE_DB.find_one({'_id': rxd_id})
			if not rxd: continue 
			if complete_only and 'complete' not in rxd:
				continue
			done_ids.add(rxd_id)
			yield (rx, rxd)

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing; defaults to False')
	parser.add_argument('-n', '--num', type = int, default = 50,
						help = 'Maximum number of records to examine; defaults to 50')
	parser.add_argument('--template_collection', type = str, default = 'transforms_forward_v1',
						help = 'Collection of templates to use; defaults to transforms_forward_v1')
	parser.add_argument('--candidate_collection', type = str, default = 'reaxys_edits_v2', 
						help = 'Collection of candidates to write to; defaults to reaxys_edits_v2')
	parser.add_argument('--mincount', type = int, default = 25,
						help = 'Minimum template count to include in transforms; defaults to 25')
	parser.add_argument('--workers', type = int, default = 10,
						help = 'Number of parallel workers, default 10')
	args = parser.parse_args()
	v = bool(args.v)

	template_collection = args.template_collection
	candidate_collection = args.candidate_collection
	n_max = int(args.num), 
	mincount = int(args.mincount)
	complete_only = True

	from rdkit import RDLogger
	lg = RDLogger.logger()
	lg.setLevel(4)

	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaxys']
	templates = db[template_collection]
	REACTION_DB = db['reactions']
	INSTANCE_DB = db['instances']
	CHEMICAL_DB = db['chemicals']
	Transformer = transformer.Transformer()
	Transformer.load(templates, mincount = mincount, get_retro = False, get_synth = True)
	print('Out of {} database templates,'.format(templates.count()))
	print('Loaded {} templates'.format(Transformer.num_templates))
	Transformer.reorder()
	print('Sorted by count, descending')
	
	db = client['prediction']
	CANDIDATE_DB = db[candidate_collection]

	# Run
	try:
		queue = Queue() 
		processes = []
		lock = Lock()
		lock_pymongo = Lock()
		for i in range(int(args.workers)):
			processes.append(
				Process(target=process_forever, args=(queue, lock_pymongo))
			)
		[p.start() for p in processes]
		print('Created {} processes'.format(len(processes)))

		# Keep queue loaded
		generator = rxd_generator()
		while True:
			if i == n_max: break
			if queue.qsize() < 100:
				queue.put(generator.next())
				i += 1
			else:
				time.sleep(1)
				continue # wait

		# Wait for queue to empty before calling .join()
		while not queue.empty():
			print('Stopped adding to queue, just waiting to empty it out...')
			time.sleep(3)

		[p.join() for p in processes]

	except KeyboardInterrupt:
		print('Stopped early, leaving pool')
