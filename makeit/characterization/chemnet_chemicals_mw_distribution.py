# Import relevant packages
from makeit.utils.chemnet_connect import *      # mongodb connection, gets 'chemicals', 'reactions'
import numpy as np     	      	   # simple calculations
import matplotlib.pyplot as plt    # visualization
import rdkit.Chem as Chem          # molecule building
import rdkit.Chem.Descriptors as Descriptors # to get MW
from makeit.utils.OPSIN_Lookup import OPSIN_Lookup # for local name-based look-up
import os                          # for saving
# Online pacckages
import cirpy					   # for resolving names/etc.
from cirpy import Molecule         # for kekulizing (better than RDKit?)

# Start OPSIN look-up
Looker = OPSIN_Lookup()

# Get distribution of (median) yields for all reactions
total_documents = chemicals.count()
all_mws = np.nan * np.ones(total_documents)
for i, chemical in enumerate(chemicals.find()):

	# Read MW
	try: 
		all_mws[i] = chemical['mol_weight']

	# No MW entry for some reason
	except: 
		all_mws[i] = np.nan

	# If no entry, need to resolve using CIRPy (tries OPSIN, PubChem)
	if np.isnan(all_mws[i]):
		smiles = ''
		done = False

		# Try with SMILES as-is
		if chemical['smiles'] and not done:
			try:
				this_mol = Chem.MolFromSmiles(chemical['smiles'])
				all_mws[i] = Descriptors.MolWt(this_mol)
				done = True
			except:
				pass

		# # Try using OPSIN to parse name
		# if chemical['name'] and not done:
		# 	try:
		# 		smiles = Looker.lookup(chemical['name'][0])
		# 		if smiles:
		# 			this_mol = Chem.MolFromSmiles(smiles)
		# 			all_mws[i] = Descriptors.MolWt(this_mol)
		# 			done = True
		# 	except:
		# 		pass

		# # ONLINE: Try with SMILES read by cirpy
		# if chemical['smiles'] and not done:
		# 	try:
		# 		smiles = cirpy.resolve(chemical['smiles'], 'smiles')
		# 		if smiles:
		# 			this_mol = Chem.MolFromSmiles(smiles)
		# 			all_mws[i] = Descriptors.MolWt(this_mol)
		# 			done = True
		# 	except:
		# 		pass

		# # ONLINE: Try with CAS
		# if chemical['cas'] and not done:
		# 	try:
		# 		smiles = cirpy.resolve(chemical['cas'], 'smiles')
		# 		if smiles:
		# 			this_mol = Chem.MolFromSmiles(smiles)
		# 			all_mws[i] = Descriptors.MolWt(this_mol)
		# 			done = True
		# 	except:
		# 		pass

		# # ONLINE: Try with name
		# if chemical['name'] and not done:
		# 	try:
		# 		smiles = cirpy.resolve(chemical['name'][0], 'smiles')
		# 		if smiles:
		# 			this_mol = Chem.MolFromSmiles(smiles)
		# 			all_mws[i] = Descriptors.MolWt(this_mol)
		# 			done = True
		# 	except:
		# 		pass

		# Report unprocessed documents
		if not done:
			print 'failed to resolve ' + str(chemical['_id'])

	if (i % 10000) == 0:
		print 'Completed %i/%i' % (i, total_documents)

	# # Terminate early (for testing)
	# if i > 312000: 
	# 	all_mws = all_mws[:312000]
	# 	break

# Filter out NaN values (which we get if we take np.median of empty array)
print 'filtering ' + str(len(all_mws[np.isnan(all_mws)])) + ' entries with mol_weight == NaN'
all_mws = all_mws[~np.isnan(all_mws)]
# Filter out large values (need to look into still...)
print 'filtering ' + str(len(all_mws[all_mws > 1500])) + ' entries with mol_weight > 1500 g/mol'
all_mws = all_mws[all_mws <= 1500]
print str(len(all_mws)) + ' entries remaining in final dataset'

# Visualize yields in histogram
weights = np.ones_like(all_mws) / len(all_mws)
n, bins, patches = plt.hist(all_mws, 50, facecolor = 'blue', alpha = 0.5, weights = weights)
plt.xlabel('Molecular weight (g/mol)')
plt.ylabel('Normalized frequency')
plt.title('Histogram of MWs for chemicals in Bishop\'s 1M reaction dataset')
plt.grid(True)
plt.savefig(os.path.dirname(__file__) + '/Histogram of chemicals mol_weight.png', bbox_inches = 'tight')