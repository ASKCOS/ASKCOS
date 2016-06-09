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

def main(template_collection = 'lowe_refs_general_v3', USE_REACTIONSMILES = True,
		mincount = 4, n_max = 50, seed = None, outfile = '.'):

	from rdkit import RDLogger
	lg = RDLogger.logger()
	lg.setLevel(4)

	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['askcos_transforms']
	templates = db[template_collection]
	Transformer = transformer.Transformer()
	Transformer.load(templates, mincount = mincount, get_retro = False, get_synth = True)
	print('Out of {} database templates,'.format(templates.count()))
	print('Loaded {} templates'.format(Transformer.num_templates))
	Transformer.reorder()
	print('Sorted by count, descending')

	# Define generator
	class Randomizer():
		def __init__(self):
			self.done_ids = []
			if seed == None:
				seed = np.random.randint(10000)
			np.random.seed(seed)
			if outfile:
				with open(os.path.join(outfile, 'seed.txt'), 'w') as fid:
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

	smilesfixer = SmilesFixer()
	try:
		for i, reaction in generator:

			if i == n_max: 
				break

			print('#########')
			print('## RXN {}'.format(i))
			print('#########')

			if bool(USE_REACTIONSMILES):
				rxn_smiles = reaction['reaction_smiles'].split(' ')[0]
				all_smiles = [smilesfixer.fix_smiles(x) for x in rxn_smiles.split('>')[0].split('.')]
				mol = Chem.MolFromSmiles(rxn_smiles.split('>')[2])
			else:
				all_smiles =  [smilesfixer.fix_smiles(x['smiles']) for x in reaction['reactants']]
				if 'catalysts' in reaction:
					all_smiles += [smilesfixer.fix_smiles(x['smiles']) for x in reaction['catalysts']] 
				if 'spectators' in reaction:
					all_smiles += [smilesfixer.fix_smiles(x['smiles']) for x in reaction['spectators']] 
				mol = Chem.MolFromSmiles(reaction['products'][0]['smiles'])
			
			if mol:
				[x.ClearProp('molAtomMapNumber') for x in mol.GetAtoms()] # remove atom mapping
				print('REACTANTS: {}'.format('.'.join(all_smiles)))
				print('PRODUCT: {}'.format(Chem.MolToSmiles(mol, isomericSmiles = USE_STEREOCHEMISTRY)))
				target_smiles = smilesfixer.fix_smiles(Chem.MolToSmiles(mol, isomericSmiles = USE_STEREOCHEMISTRY))
				result = Transformer.perform_forward('.'.join(all_smiles), progbar = True)

			found_true = False
			for product in result.products:
				if target_smiles in product.smiles_list:
					found_true = True
			if found_true:
				print('True product found!')
			else:
				print('True product not found...')
				raw_input('Pausing...')
	except Exception as e:
		print('Error! {}'.format(e))




if __name__ == '__main__':
	main()