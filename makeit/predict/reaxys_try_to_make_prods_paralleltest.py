# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import argparse
import numpy as np     	      	   # for simple calculations
import matplotlib
import matplotlib.pyplot as plt    # for visualization
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import os                          # for saving
matplotlib.rc('font', **{'size': 18})
from makeit.retro.draw import ReactionStringToImage
from makeit.retro.canonicalization import SmilesFixer
import sys
import re
import time
from tqdm import tqdm
from multiprocessing import Pool
import itertools
from functools import partial

def apply_this_template(this_template, reactants):
	'''Given an RDKit reaction and a reactants molecule, apply
	the template and return a set of outcome SMILES strings'''
	prods = set()
	# Perform
	try:
		outcomes = this_template.RunReactants([reactants])
	except Exception as e:
		#print(e)
		return []
	if not outcomes: return []
	for j, outcome in enumerate(outcomes):
		try:
			for x in outcome:
				#print('  - {}'.format(Chem.MolToSmiles(x)))
				x.UpdatePropertyCache()
				Chem.SanitizeMol(x)
		except Exception as e:
			#print(e)
			continue
		product_smiles = [Chem.MolToSmiles(x, isomericSmiles = USE_STEREOCHEMISTRY) for x in outcome]
		prods.add('.'.join(sorted(product_smiles)))
	return prods 

def multi_arg_apply(args):
	'''Because Pool must use a single arg function defined at the highest
	level (importable from __main__), we must use this wrapper'''
	return apply_this_template(*args)

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing (incl. saving images); defaults to False')
	parser.add_argument('-o', '--out', type = str, default = 'test/rxn_test',
						help = 'Folder to output images to if defined')
	parser.add_argument('-n', '--num', type = int, default = 50,
						help = 'Maximum number of records to examine; defaults to 50')
	parser.add_argument('-t', '--templates', type = str, default = 'transforms_forward_v1',
						help = 'Name of reaction template collection, defaults to transforms_forward_v1')
	parser.add_argument('--mincount', type = int, default = 10, 
						help = 'Minimum count of reaction templates to use, defaults to 2')
	parser.add_argument('--seed', type = int, default = None,
						help = 'Seed for instance randomization')
	parser.add_argument('--skip_remaining', type = str, default = 'n',
						help = 'Whether to skip remaining templates when true found, default n')
	parser.add_argument('-w', '--workers', type = int, default = 2,
						help = 'Number of workers in parallel pool, default 2')
	parser.add_argument('-c', '--chunksize', type = int, default = 1,
						help = 'Chunk size for parallel pool, default 1')

	args = parser.parse_args()

	v = args.v
	from rdkit import RDLogger
	lg = RDLogger.logger()
	if not v: 
		lg.setLevel(RDLogger.CRITICAL)
	from rdkit import rdBase
	rdBase.DisableLog('rdApp.error')

	N = int(args.num)
	skip_remaining = args.skip_remaining in ['y', 'Y', 'yes', 'Yes', '1', 't', 'True', 'T']
	N_THREADS = int(args.workers)
	chunksize = int(args.chunksize)

	# DATABASE
	from pymongo import MongoClient    # mongodb plugin
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaxys']
	reactions = db['reactions']
	instances = db['instances']
	chemicals = db['chemicals']

	templates = db[args.templates]
	import makeit.retro.transformer as transformer 
	docs = sorted(templates.find({'count': {'$gte': args.mincount}}, ['count', 'reaction_smarts']), key = lambda x: x['count'])
	all_reaction_smarts = [doc['reaction_smarts'] for doc in docs] 
	template_rxns = []
	for reaction_smarts in all_reaction_smarts:
		try:
			# Load as a unimolecular reaction
			unimol = '(' + reaction_smarts.split('>')[0] + ')>>(' + reaction_smarts.split('>')[2] + ')'
			rxn = AllChem.ReactionFromSmarts(str(unimol))
			if rxn.Validate()[1] == 0:
				template_rxns.append(rxn)
			else:
				print('Error loading template: {}'.format(reaction_smarts))
		except Exception as e:
			print(e)
	print('Out of {} database templates,'.format(templates.count()))
	print('Loaded {} templates'.format(len(template_rxns)))
	print('Sorted by count, descending')


	# Define generator to grab XRDs
	class Randomizer():
		def __init__(self):
			self.done_ids = []
			if args.seed:
				seed = int(args.seed)
			else:
				seed = np.random.randint(10000)
			np.random.seed(seed)
		def get_rand(self):
			'''Random WITHOUT replacement'''
			# while True:
			# 	rx = reactions.find({'random': { '$gte': np.random.random()}}).sort('random', 1).limit(1)[0]
			for rx in reactions.find():
				if rx['RX_PXRN'] and rx['RX_RXRN']: # if reactants and products recorded
					for rxd in instances.find({'RX_ID': rx['_id']}):
						if rxd['_id'] in self.done_ids: continue 
						yield (rx, rxd)
						self.done_ids.append(rxd['_id'])
	randomizer = Randomizer()
	generator = enumerate(randomizer.get_rand())

	rxn_successful = 0; rxn_unsuccessful = 0;
	smilesfixer = SmilesFixer()
	cum_time = 0
	p = Pool(N_THREADS)
	try:
		for i, (rx, rxd) in generator:

			start_time = time.time()

			if i == N: 
				N = i
				break

			print('#########')
			print('## RXN {}'.format(rxd['_id']))
			print('#########')

			rxn_smiles = rx['RXN_SMILES']
			reactants = rxn_smiles.split('>')[0]
			print('"REACTANTS": {}'.format(reactants))
			products = rxn_smiles.split('>')[2]

			# Add reagents/solvent/catalyst into reactants
			XRNs = rxd['RXD_RGTXRN'] + rxd['RXD_SOLXRN'] + rxd['RXD_CATXRN']
			if XRNs:
				docs = [chemicals.find_one({'_id': xrn}) for xrn in XRNs]
				extra_smiles = [doc['SMILES'] for doc in docs if doc]
				reactants += '.' + '.'.join([extra_smile for extra_smile in extra_smiles if extra_smile])

			reactants = Chem.MolFromSmiles(reactants)
			products = Chem.MolFromSmiles(products)

			if (not reactants) or (not products):
				print('Could not load reactants/products?')
				print('ID: {}'.format(rxd['_id']))
				continue

			[a.ClearProp('molAtomMapNumber') for a in reactants.GetAtoms()] # remove atom mapping
			[a.ClearProp('molAtomMapNumber') for a in products.GetAtoms()] # remove atom mapping

			major_product_smiles = max(Chem.MolToSmiles(products, isomericSmiles = USE_STEREOCHEMISTRY).split('.'), key = len)
			if major_product_smiles in Chem.MolToSmiles(reactants, isomericSmiles = USE_STEREOCHEMISTRY).split('.'):
				print('WARNING: recorded product already in reactants?')
				if v: raw_input('Pause...enter will skip')
				continue

			print('REACTANTS: {}'.format(Chem.MolToSmiles(reactants, isomericSmiles = USE_STEREOCHEMISTRY)))
			print('PRODUCTS:  {}'.format(Chem.MolToSmiles(products, isomericSmiles = USE_STEREOCHEMISTRY)))
			if v: raw_input('continue to try to match this rxn...')

			# Apply all templates
			unique_results = set(itertools.chain.from_iterable(
				p.imap_unordered(
					multi_arg_apply,
					itertools.izip(template_rxns, itertools.repeat(reactants)),
					chunksize = chunksize
				)
			))
			try: 
				unique_results.remove(None)
			except KeyError:
				pass

			found = any([major_product_smiles in prod for prod in unique_results])

			if found:
				rxn_successful += 1
				if not skip_remaining: print('Found true product among {} unique products'.format(len(unique_results)))
			else:
				rxn_unsuccessful += 1
				print('Did not find true product!')
				print('But found {} unique other products'.format(len(unique_results)))
				print('REACTANTS: {}'.format(Chem.MolToSmiles(reactants, isomericSmiles = USE_STEREOCHEMISTRY)))
				print('PRODUCTS:  {}'.format(Chem.MolToSmiles(products, isomericSmiles = USE_STEREOCHEMISTRY)))

				# Also save
				if args.out:
					img = ReactionStringToImage(rxn_smiles)
					img.save(os.path.join(args.out, '{}_failed_{}.png'.format(i, str(rxd['_id']))))

			print('TOTALS:')
			print('  {} successful'.format(rxn_successful))
			print('  {} unsuccessful'.format(rxn_unsuccessful))

			this_time = time.time() - start_time
			print('This example processed in {} seconds'.format(this_time))
			cum_time += this_time

	except KeyboardInterrupt:
		print('terminated early')
		N = i

	del p

	print('Out of {} reactions:'.format(N))
	print('  {} successful'.format(rxn_successful))
	print('  {} unsuccessful'.format(rxn_unsuccessful))
	print('Took {} s to calculate using {} workers'.format(cum_time, N_THREADS))
