# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import argparse
import numpy as np     	      	   # for simple calculations
import matplotlib
import matplotlib.pyplot as plt    # for visualization
import rdkit.Chem as Chem
import os                          # for saving
matplotlib.rc('font', **{'size': 18})
from makeit.retro.draw import ReactionStringToImage
from makeit.retro.canonicalization import SmilesFixer
import sys
import re

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing (incl. saving images); defaults to False')
	parser.add_argument('-o', '--out', type = str, default = 'test/rxn_test',
						help = 'Folder to output images to if defined')
	parser.add_argument('-n', '--num', type = int, default = 50,
						help = 'Maximum number of records to examine; defaults to 50')
	parser.add_argument('-e', '--examples', type = str, default = 'lowe_1976-2013_USPTOgrants_reactions',
						help = 'Name of reaction example collection, defaults to Lowe USPTO grants reactions')
	parser.add_argument('-t', '--templates', type = str, default = 'lowe_refs_general',
						help = 'Name of reaction template collection, defaults to lowe_refs_general')
	parser.add_argument('--mincount', type = int, default = 2, 
						help = 'Minimum count of reaction templates to use, defaults to 2')
	parser.add_argument('-s', '--shuffle', type = bool, default = True,
						help = 'Shuffle reaction example database? defaults to True')
	parser.add_argument('--seed', type = int, default = None,
						help = 'Seed for random number generator')

	args = parser.parse_args()

	v = args.v
	from rdkit import RDLogger
	lg = RDLogger.logger()
	if not v: lg.setLevel(4)

	N = int(args.num)

	# DATABASE
	from pymongo import MongoClient    # mongodb plugin
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaction_examples']
	reactions = db[args.examples]

	import makeit.retro.transformer as transformer 
	db = client['askcos_transforms']
	templates = db[args.templates]
	Transformer = transformer.Transformer()
	Transformer.load(templates, mincount = args.mincount, get_retro = False, get_synth = True)
	print('Out of {} database templates,'.format(templates.count()))
	print('Loaded {} templates'.format(Transformer.num_templates))
	Transformer.reorder()
	print('Sorted by count, descending')

	# Iterate through all examples and see if we made it


	# Define generator
	if bool(args.shuffle):
		class Randomizer():
			def __init__(self):
				self.done_ids = []
				if args.seed:
					seed = int(args.seed)
				else:
					seed = np.random.randint(10000)
				np.random.seed(seed)
				if args.out:
					with open(os.path.join(args.out, 'seed.txt'), 'w') as fid:
						fid.write('{}'.format(seed))
			def get_rand(self):
				'''Random WITHOUT replacement'''
				while True:
					doc = reactions.find({'products': {'$size': 1}, \
						'products.0.smiles': {'$not': re.compile('.*\..*')},
						'random': { '$gte': np.random.random()}}).sort('random', 1).limit(1)[0]
					if doc['_id'] in self.done_ids: continue
					self.done_ids.append(doc['_id'])
					yield doc
		randomizer = Randomizer()
		generator = enumerate(randomizer.get_rand())
	else:
		generator = enumerate(reactions.find({'products': {'$size': 1}, \
			'products.0.smiles': {'$not': re.compile('.*\..*')}}, no_cursor_timeout=True))


	rxn_successful = 0; rxn_unsuccessful = 0;
	smilesfixer = SmilesFixer()
	try:
		for i, reaction in generator:
			if i == N: 
				N = i
				break

			print('#########')
			print('## RXN {}'.format(i))
			print('#########')

			all_smiles =  [smilesfixer.fix_smiles(x['smiles']) for x in reaction['reactants']]
			if 'catalysts' in reaction:
				all_smiles += [smilesfixer.fix_smiles(x['smiles']) for x in reaction['catalysts']] 
			if 'spectators' in reaction:
				all_smiles += [smilesfixer.fix_smiles(x['smiles']) for x in reaction['spectators']] 

			mol = Chem.MolFromSmiles(reaction['products'][0]['smiles'])
			if mol:
				product_smiles = Chem.MolToSmiles(mol, isomericSmiles = USE_STEREOCHEMISTRY)
				result = Transformer.perform_forward('.'.join(all_smiles), 
					stop_if = smilesfixer.fix_smiles(product_smiles))
			else: 
				result = False

			if result == True:
				rxn_successful += 1
			else:
				rxn_unsuccessful += 1
				# Also save
				print(reaction)
				if args.out:
					img = ReactionStringToImage(reaction['reaction_smiles'].split(' ')[0])
					img.save(os.path.join(args.out, '{}_failed_{}.png'.format(i, str(reaction['_id']))))

	except KeyboardInterrupt:
		print('terminated early')
		N = i

	print('Out of {} reactions:'.format(N))
	print('  {} successful'.format(rxn_successful))
	print('  {} unsuccessful'.format(rxn_unsuccessful))