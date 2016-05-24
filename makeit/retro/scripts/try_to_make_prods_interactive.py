# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import argparse
import numpy as np     	      	   # for simple calculations
import matplotlib
import matplotlib.pyplot as plt    # for visualization
import os                          # for saving
matplotlib.rc('font', **{'size': 18})
from makeit.retro.draw import ReactionStringToImage, MolsSmilesToImage
import sys

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing (incl. saving images); defaults to False')
	parser.add_argument('-o', '--out', type = str, default = 'test/rxn_test',
						help = 'Folder to output images to if defined')
	parser.add_argument('-t', '--templates', type = str, default = 'lowe_refs_general_v2',
						help = 'Name of reaction template collection, defaults to lowe_refs_general_v2')
	parser.add_argument('--mincount', type = int, default = 2, 
						help = 'Minimum count of reaction templates to use, defaults to 2')

	args = parser.parse_args()

	v = args.v
	from rdkit import RDLogger
	lg = RDLogger.logger()
	if not v: lg.setLevel(4)

	import makeit.retro.transformer as transformer 
	from pymongo import MongoClient    # mongodb plugin
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['askcos_transforms']
	templates = db[args.templates]
	Transformer = transformer.Transformer()
	Transformer.load(templates, mincount = args.mincount, get_retro = False, get_synth = True)
	print('Out of {} database templates,'.format(templates.count()))
	print('Loaded {} templates'.format(Transformer.num_templates))
	Transformer.reorder()
	print('Sorted by count, descending')

	done = False
	n_top = 5
	while not done:
		reactant_smiles_string = raw_input('Enter reactant SMILES strings: ')
		if reactant_smiles_string in ['done', 'quit', 'exit']:
			done = True
			continue
		if 'n = ' in reactant_smiles_string:
			n_top = int(reactant_smiles_string.strip().split(' ')[-1])
			continue


		try:
			try:
				out_folder = os.path.join(args.out, reactant_smiles_string)
				os.mkdir(out_folder)
			except:
				pass
			result = Transformer.perform_forward(reactant_smiles_string.strip())
			img = ReactionStringToImage(reactant_smiles_string + '>>')
			img.save(os.path.join(out_folder, 'reactants.png'))

			for product in result.return_top(n = n_top):
				img = MolsSmilesToImage(product['smiles'])
				img.save(os.path.join(out_folder, 
					'rank{}_score{}.png'.format(product['rank'], product['num_examples'])))

			print('Total of {} possible products'.format(len(result.products)))

		except Exception as e:
			print(e)