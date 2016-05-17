# Import relevant packages
import numpy as np     	      	   # for simple calculations
import matplotlib
import matplotlib.pyplot as plt    # for visualization
import os                          # for saving
matplotlib.rc('font', **{'size': 18})
from makeit.retro.draw import ReactionStringToImage
import sys

# NUMBER TO TRY
N = int(sys.argv[1])

# SILENCE RDKIT WARNINGS
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(4)

# DATABASE
from pymongo import MongoClient    # mongodb plugin
client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaction_examples']
reactions = db['lowe_1976-2013_USPTOgrants_reactions']

import makeit.retro.transformer as transformer 
db = client['askcos_transforms']
templates = db['lowe_refs_general']
Transformer = transformer.Transformer()
Transformer.load(templates)
print('Loaded {} templates'.format(Transformer.num_templates))

# Iterate through all examples and see if we made it
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
			print(reaction)

except KeyboardInterrupt:
	print('terminated early')
	N = i

print('Out of {} reactions:'.format(N))
print('  {} successful'.format(rxn_successful))
print('  {} unsuccessful'.format(rxn_unsuccessful))