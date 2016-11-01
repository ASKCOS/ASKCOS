# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import argparse
import numpy as np     	      	   # for simple calculations
import os                          # for saving
import sys
import rdkit.Chem as Chem
import makeit.retro.transformer as transformer 
from makeit.retro.canonicalization import SmilesFixer
from pymongo import MongoClient    # mongodb plugin
from makeit.predict.summarize_reaction_outcome import summarize_reaction_outcome
from makeit.embedding.descriptors import edits_to_vectors, oneHotVector # for testing
import re
import time
from tqdm import tqdm
from multiprocessing import Process, Lock, Queue
# Define generator
class Randomizer():
	def __init__(self, seed, done_ids = []):
		self.done_ids = done_ids
		np.random.seed(seed)
		with open(os.path.join(os.path.dirname(__file__), 'reaxys_generate_candidate_edits_seed.txt'), 'w') as fid:
			fid.write('{}'.format(seed))

	def get_sequential(self):
		for rx in REACTION_DB.find({'_id': {'$gt': start_at_id}}, no_cursor_timeout = True):
			if not rx: break 
			for i in range(1, rx['RX_NVAR'] + 1):
				rxd_id = '{}-{}'.format(rx['_id'], i)
				if rxd_id in self.done_ids: continue
				rxd = INSTANCE_DB.find_one({'_id': rxd_id})
				if not rxd: continue 
				if complete_only and 'complete' not in rxd:
					continue
				self.done_ids.append(rxd_id)
				yield (rx, rxd)


def process_one(queue, lock_pymongo):
	i = 0
	(rx, rxd) = queue.get()

	try:

		# LOGGING
		start_time = time.time()

		if i == n_max: 
			raise KeyboardInterrupt

		print('#########')
		print('## RXN {}'.format(rxd['_id']))
		print('#########')

		rxn_smiles = rx['RXN_SMILES']
		reactants = rxn_smiles.split('>')[0]
		products = rxn_smiles.split('>')[2]

		# Add reagents/solvent/catalyst into reactants
		XRNs = rxd['RXD_RGTXRN'] + rxd['RXD_SOLXRN'] + rxd['RXD_CATXRN']
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
			if v: raw_input('Pause...enter will skip')
			return

		#print('REACTANTS: {}'.format(Chem.MolToSmiles(reactants, isomericSmiles = USE_STEREOCHEMISTRY)))
		#print('PRODUCTS:  {}'.format(Chem.MolToSmiles(products, isomericSmiles = USE_STEREOCHEMISTRY)))
		if v: raw_input('continue to try to match this rxn...')


		n_reactant_atoms = len(reactants.GetAtoms())
		#print('Number of reactant atoms: {}'.format(n_reactant_atoms))
		if n_reactant_atoms > 80:
			#print('Skipping huge molecule! N_reactant_atoms = {}'.format(n_reactant_atoms))
			return

		# Report current reactant SMILES string
		#print('Reactants w/o map: {}'.format(Chem.MolToSmiles(reactants)))
		# Add new atom map numbers
		[a.SetProp('molAtomMapNumber', str(i+1)) for (i, a) in enumerate(reactants.GetAtoms())]
		# Report new reactant SMILES string
		#print('Reactants w/ map: {}'.format(Chem.MolToSmiles(reactants)))

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
					[a.SetProp('molAtomMapNumber', a.GetProp('old_molAtomMapNumber')) \
						for (i, a) in enumerate(outcome.GetAtoms()) \
						if 'old_molAtomMapNumber' in a.GetPropsAsDict()]
				except Exception as e:
					if v: print(e)
					continue
				if v: print('Outcome SMILES: {}'.format(Chem.MolToSmiles(outcome)))

				# Reduce to largest (longest) product only?
				candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY)
				if singleonly: 
					candidate_smiles = max(candidate_smiles.split('.'), key = len)
					outcome = Chem.MolFromSmiles(candidate_smiles)
					
				# Find what edits were made
				edits = summarize_reaction_outcome(reactants, outcome)
				if v: print(edits)

				# Remove mapping before matching
				[x.ClearProp('molAtomMapNumber') for x in outcome.GetAtoms() if x.HasProp('molAtomMapNumber')] # remove atom mapping from outcome
				if Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY) == major_product_smiles:
					if v: print('Matched true [{}]'.format(major_product_smiles))
					found_true = True

				# Overwrite candidate_smiles without atom mapping numbers
				candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY)

				# Add to ongoing list
				if (candidate_smiles, edits) not in candidate_edits:
					candidate_edits.append((candidate_smiles, edits))

		# Prepare doc and insert
		if found_true: print('Found true product for {}'.format(rxd['_id']))
		doc = {
			'_id': rxd['_id'],
			'reactant_smiles': Chem.MolToSmiles(reactants, isomericSmiles = USE_STEREOCHEMISTRY),
			'edit_candidates': candidate_edits,
			'product_smiles_true': major_product_smiles,
			'found': found_true,
			'num_candidates': len(candidate_edits),
		}
		try:
			lock_pymongo.acquire()
			res = candidates.insert(doc)
			lock_pymongo.release()
		except Exception as e:
			print(e)
			return

		# LOGGING
		end_time = time.time()
		print('time: {}'.format(end_time - start_time))
		print('unique edit sets using longest prod: {}'.format(len(candidate_edits)))
		if log: flog.write('{}\t{}\t{}\t{}\n'.format(i, n_reactant_atoms, len(candidate_edits), end_time - start_time))

	except KeyboardInterrupt:
		print('Breaking early')
		raise KeyboardInterrupt
	# except Exception as e:
	# 	print('Error, {}'.format(e))
	# 	return


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing; defaults to False')
	parser.add_argument('-n', '--num', type = int, default = 50,
						help = 'Maximum number of records to examine; defaults to 50')
	parser.add_argument('--template_collection', type = str, default = 'transforms_forward_v1',
						help = 'Collection of templates to use; defaults to transforms_forward_v1')
	parser.add_argument('--candidate_collection', type = str, default = 'reaxys_edits_v1', 
						help = 'Collection of candidates to write to; defaults to reaxys_edits_v1')
	parser.add_argument('--seed', type = int, default = None,
						help = 'Seed for random number generator')
	parser.add_argument('--mincount', type = int, default = 100,
						help = 'Minimum template count to include in transforms; defaults to 100')
	parser.add_argument('--singleonly', type = bool, default = True,
						help = 'Whether to record major product only; defaults to True')
	parser.add_argument('--check', type = bool, default = True,
						help = 'Whether to check current collection to see if reaction example has been done')
	parser.add_argument('--log', type = bool, default = True,
						help = 'Whether to log wall times / number of candidate atoms / etc., default True')
	parser.add_argument('--complete_only', type = str, default = 'y',
						help = 'Whether to only use complete reaction instances, default y')
	parser.add_argument('--workers', type = int, default = 10,
						help = 'Number of parallel workers, default 10')
	args = parser.parse_args()
	v = bool(args.v)

	template_collection = args.template_collection
	candidate_collection = args.candidate_collection
	seed = args.seed
	n_max = int(args.num), 
	mincount = int(args.mincount)
	singleonly = bool(args.singleonly),
	check = bool(args.check)
	log = bool(args.log)
	complete_only = args.complete_only in ['y', 'Y', 'true', 'T', '1']
	if complete_only:
		print('Using only complete reaction instances')

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
	candidates = db[candidate_collection]

	if check:
		done_ids = [doc['_id'] for doc in candidates.find(
			{}, ['_id']
		)]
		print('Checked completed entries')

	else:
		done_ids = []

	start_at_id = max([int(rxd_id.split('-')[0]) for rxd_id in done_ids]) - 1
	print('Starting at RX_ID {}'.format(start_at_id))

	if seed == None:
		seed = np.random.randint(100000)
	else:
		seed = int(seed)
	randomizer = Randomizer(seed, done_ids = done_ids)
	generator = enumerate(randomizer.get_sequential())

	smilesfixer = SmilesFixer()

	# LOGGING
	if log:
		flog = open('GENERATE_CANDIDATES_LOG_{}.txt'.format(seed), 'w')
		flog.write('mincount: {}\n'.format(mincount))
		flog.write('number of templates: {}\n'.format(Transformer.num_templates))

	# for i, (rx, rxd) in generator:


	# Run
	try:
		queue = Queue() 
		processes = []
		lock = Lock()
		lock_pymongo = Lock()
		for i in range(int(args.workers)):
			processes.append(
				Process(target=process_one, args=(queue, lock_pymongo))
			)
		[p.start() for p in processes]

		# Add to queue now
		for i, (rx, rxd) in generator:
			if i == n_max: break
			while True:
				if queue.qsize() < 100:
					queue.put((rx, rxd))
					break
				else:
					continue # wait
		print('Added {} examples to queue'.format(i+1))
		#

		[p.join() for p in processes]
	except KeyboardInterrupt:
		print('Stopped early, leaving pool')

	if log: flog.close()