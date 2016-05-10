# Import relevant packages
import numpy as np     	      	   # for simple calculations
import matplotlib
import matplotlib.pyplot as plt    # for visualization
import os                          # for saving
matplotlib.rc('font', **{'size': 18})

# DATABASE
from pymongo import MongoClient    # mongodb plugin
client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaction_examples']
reactions = db['lowe_1976-2013_USPTOgrants_reactions']

import makeit.retro.transformer as transformer 
db = client['askcos_transforms']
templates = db['lowe']
Transformer = transformer.Transformer()
Transformer.load(templates)
print('Loaded {} templates'.format(Transformer.num_templates))

# Iterate through all examples and see if we made it
N = reactions.count()
rxn_successful = 0; rxn_unsuccessful = 0;
for i, reaction in enumerate(reactions.find({'products': {'$size': 1}})):
	if i == 10: 
		N = i
		break

	all_smiles =  [x['smiles'] for x in reaction['reactants']]
	if 'catalysts' in reaction:
		all_smiles += [x['smiles'] for x in reaction['catalysts']] 
	if 'spectators' in reaction:
		all_smiles += [x['smiles'] for x in reaction['spectators']] 

	result = Transformer.perform_forward('.'.join(all_smiles))

	success = False
	for product in sorted(result.products, key = lambda x: x.num_examples, reverse = True):
		if '.'.join(product.smiles_list) == reaction['products'][0]['smiles']:
			success = True
			break
	if success:
		rxn_successful += 1
	else:
		rxn_unsuccessful += 1
		print reaction

print('Out of {} reactions:'.format(N))
print('  {} successful'.format(rxn_successful))
print('  {} unsuccessful'.format(rxn_unsuccessful))