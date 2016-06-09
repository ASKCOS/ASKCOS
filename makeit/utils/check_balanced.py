## Meant to use with Daniel Lowe's patent database
# reaction SMILES strings only (for now)

from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import argparse
from numpy.random import shuffle # for random selection
import rdkit.Chem as Chem          # molecule building
from rdkit.Chem import AllChem
from collections import defaultdict
import rdkit.Chem.Draw as Draw
from rdkit import RDLogger
import datetime # for info files
import json # for dumping
import sys  # for commanad line
import os   # for file paths
import re 
import itertools
from makeit.retro.draw import *


def main(N = 15):
	'''Read reactions from Lowe's patent reaction SMILES'''
	global v
	# Open file

	# Define scoring variables
	total_attempted = 0 # total reactions simulated (excludes skipped)
	total_balanced = 0 # actual products predicted
	total_unbalanced = 0 # ONLY actual products predicted

	N = min([N, example_collection.count()])

	try: # to allow breaking
		# Look for entries
		i = -1
		for example_doc in example_collection.find(no_cursor_timeout=True):
			i += 1
			if i < 616120: continue
			# Are we done?
			if i == N:
				i -= 1
				break

			# KNOWN ISSUES
			if i > 13700 and i < 13800: continue
			if i > 23283 and i < 23285: continue
			if i > 23860 and i < 23863: continue
			if i > 56778 and i < 56782: continue
			if i > 87600 and i < 87700: continue
			if i > 88200 and i < 88300: continue
			if i > 185000 and i < 185100: continue

			if v: print('### Reaction {}'.format(i))

			try:
				# Unpack
				reaction_smiles = str(example_doc['reaction_smiles'])
				reactants, agents, products = [mols_from_smiles_list(x) for x in 
											[mols.split('.') for mols in reaction_smiles.split('>')]]
				[Chem.SanitizeMol(mol) for mol in reactants + agents + products]
			except:
				# can't sanitize -> skip
				if v: print('- skipping -')
				continue

			all_reactant_symbols = []
			all_product_symbols = []
			reactant_Hs = 0
			product_Hs = 0
			for reactant in reactants:
				all_reactant_symbols += [a.GetSymbol() for a in reactant.GetAtoms()]
				reactant_Hs += sum([a.GetTotalNumHs() for a in reactant.GetAtoms()])
			for product in products:
				all_product_symbols += [a.GetSymbol() for a in product.GetAtoms()]
				product_Hs += sum([a.GetTotalNumHs() for a in product.GetAtoms()])

			balanced = True
			if reactant_Hs != product_Hs:
				if v: print('- different number of hydrogens! {} v. {}'.format(reactant_Hs, product_Hs))
				balanced = False 
			if len(all_reactant_symbols) != len(all_product_symbols):
				if v: print('- different number of heavy atoms! {} v. {}'.format(len(all_reactant_symbols), len(all_product_symbols)))
				balanced = False
			elif sorted(all_reactant_symbols) != sorted(all_product_symbols):
				if v: print('- different atomic identitites (unexpected...)')
				balanced = False 

			if balanced:
				total_balanced += 1
				example_collection.update_one(
					{'_id': example_doc['_id']},
					{
						'$set': {
							'balanced': True,
						},
					}
				)
			else:
				total_unbalanced += 1
			total_attempted += 1

			del all_reactant_symbols
			del all_product_symbols

			# Report progress
			if (i % 1000) == 0:
				print('{}/{}'.format(i, N))
				total_examples = i + 1
				print('...finished looking through {} reaction records'.format(min([N, i + 1])))

				print('Error-free parsing in {}/{} ({}%) cases'.format(total_attempted, 
					total_examples, total_attempted * 100.0 / total_examples))
				print('Balanced reactions in {}/{} ({}%) cases'.format(total_balanced, 
					total_examples, total_balanced * 100.0 / total_examples))
				print('Unbalanced reactions in {}/{} ({}%) cases'.format(total_unbalanced, 
					total_examples, total_unbalanced * 100.0 / total_examples))


			# Pause
			#if v: raw_input('Enter anything to continue...')

	except KeyboardInterrupt:
		print('Stopped early!')		
	except Exception as e:
		print(e)

	total_examples = i + 1
	print('...finished looking through {} reaction records'.format(min([N, i + 1])))

	print('Error-free parsing in {}/{} ({}%) cases'.format(total_attempted, 
		total_examples, total_attempted * 100.0 / total_examples))
	print('Balanced reactions in {}/{} ({}%) cases'.format(total_balanced, 
		total_examples, total_balanced * 100.0 / total_examples))
	print('Unbalanced reactions in {}/{} ({}%) cases'.format(total_unbalanced, 
		total_examples, total_unbalanced * 100.0 / total_examples))
	

	return True


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing; defaults to False')
	parser.add_argument('-n', '--num', type = int, default = 50,
						help = 'Maximum number of records to examine; defaults to 50')
	parser.add_argument('-c', '--collection', type = str, default = 'lowe_1976-2013_USPTOgrants',
						help = 'Collection in reaction_examples to use; defaults to lowe_1976-2013_USPTOgrants')
	args = parser.parse_args()


	# DATABASE
	from pymongo import MongoClient    # mongodb plugin
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaction_examples']
	example_collection = db[args.collection]

	v = args.v
	lg = RDLogger.logger()
	if not v: lg.setLevel(4)

	main(N = args.num)