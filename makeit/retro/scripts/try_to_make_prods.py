# Import relevant packages
from __future__ import print_function
import argparse
import numpy as np     	      	   # for simple calculations
import matplotlib
import matplotlib.pyplot as plt    # for visualization
import os                          # for saving
matplotlib.rc('font', **{'size': 18})
from makeit.retro.draw import ReactionStringToImage
import sys

# NUMBER TO TRY
N = int(sys.argv[1])

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing (incl. saving images); defaults to False')
	parser.add_argument('-o', '--out', type = str, default = 'test/rxn_test',
						help = 'Folder to output images to if defined')
	parser.add_argument('-n', '--num', type = int, default = 50,
						help = 'Maximum number of records to examine; defaults to 50')
	parser.add_argument('-e', '--examples', type = str, default = 'lowe_1976-2013_USPTOgrants_reactions',
						help = 'Name of reaction example collection, defaults to Lowe USPTO grants')
	parser.add_argument('-t', '--templates', type = str, default = 'lowe_refs_general',
						help = 'Name of reaction template collection, defaults to lowe_refs_general')
	parser.add_argument('--mincount', type = int, default = 2, 
						help = 'Minimum count of reaction templates to use, defaults to 2')
	# parser.add_argument('-s', '--shuffle', type = bool, default = True,
	# 					help = 'Shuffle reaction example database? defaults to True')

	args = parser.parse_args()

	v = args.v
	from rdkit import RDLogger
	lg = RDLogger.logger()
	if not v: lg.setLevel(4)

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

	# Iterate through all examples and see if we made it
	# matched_docs = reactions.find({'products': {'$size': 1}}, no_cursor_timeout=True)
	# ids = [x['_id'] for x in matched_docs]
	# if bool(args.shuffle):
	# 	np.shuffle(ids)	
	rxn_successful = 0; rxn_unsuccessful = 0;
	try:
		for i, reaction in enumerate(reactions.find({'products': {'$size': 1}}, no_cursor_timeout=True)):
			# if i < 64: continue
			#if i < 1000: continue
			if i == N: 
				N = i
				break

			print('#########')
			print('## RXN {}'.format(i))
			print('#########')

			all_smiles =  [x['smiles'] for x in reaction['reactants']]
			if 'catalysts' in reaction:
				all_smiles += [x['smiles'] for x in reaction['catalysts']] 
			if 'spectators' in reaction:
				all_smiles += [x['smiles'] for x in reaction['spectators']] 

			result = Transformer.perform_forward('.'.join(all_smiles), stop_if = reaction['products'][0]['smiles'])
			if result == True:
				rxn_successful += 1
			else:
				rxn_unsuccessful += 1
				# Also save
				print(reaction)
				if args.out:
					img = ReactionStringToImage(reaction['reaction_smiles'].split(' ')[0])
					img.save(os.path.join(args.out, '{}_failed.png'.format(i)))

	except KeyboardInterrupt:
		print('terminated early')
		N = i

	print('Out of {} reactions:'.format(N))
	print('  {} successful'.format(rxn_successful))
	print('  {} unsuccessful'.format(rxn_unsuccessful))