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
from tqdm import tqdm

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

	args = parser.parse_args()

	v = args.v
	from rdkit import RDLogger
	lg = RDLogger.logger()
	if not v: 
		lg.setLevel(RDLogger.CRITICAL)
	from rdkit import rdBase
	rdBase.DisableLog('rdApp.error')

	N = int(args.num)

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
	try:
		for i, (rx, rxd) in generator:

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


			found = False
			unique_results = set()
			for template_rxn in tqdm(template_rxns):
				# Perform
				try:
					outcomes = template_rxn.RunReactants([reactants])
				except Exception as e:
					#print(e)
					continue
				if not outcomes: continue
				#print('TEMPLATE: {}'.format(AllChem.ReactionToSmarts(template_rxn)))
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
					if major_product_smiles in product_smiles:
						found = True 
						break
					unique_results.add('.'.join(sorted(product_smiles)))
				if found: 
					print('Found true product, total of {} unique candidates - skipping remaining templates'.format(len(unique_results)))
					break
			if found:
				rxn_successful += 1
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

	except KeyboardInterrupt:
		print('terminated early')
		N = i

	print('Out of {} reactions:'.format(N))
	print('  {} successful'.format(rxn_successful))
	print('  {} unsuccessful'.format(rxn_unsuccessful))
