# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import argparse
import numpy as np     	      	   # for simple calculations
import os                          # for saving
import sys
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from pymongo import MongoClient    # mongodb plugin
from makeit.predict.summarize_reaction_outcome import summarize_reaction_outcome

def main(TRANSFORM_DB, v = False):

	from rdkit import RDLogger
	lg = RDLogger.logger()
	lg.setLevel(4)
	
	print('# templates in database: {}'.format(TRANSFORM_DB.count()))
	db_filter = {'edits': {'$exists': False}}
	db_filter = {}
	N_left = TRANSFORM_DB.count(db_filter)
	print('# templates in database w/o edits: {}'.format(N_left))
	for i, template in enumerate(TRANSFORM_DB.find(db_filter)):
		print('Template {}/{}'.format(i+1, N_left))

		try:
			rxn_string = template['rxn_example']
			rxn_smarts_split = str(template['reaction_smarts']).split('>') # retro
			rxn_smarts = '(' + rxn_smarts_split[2] + ')>>(' + rxn_smarts_split[0] + ')' # forward
			rxn = AllChem.ReactionFromSmarts(rxn_smarts)
			rxn.Validate()

			reactants = Chem.MolFromSmiles(rxn_string.split('>')[0])
			products = Chem.MolFromSmiles(rxn_string.split('>')[2])

			reactants = rxn.GetReactants()[0]
			reactants.UpdatePropertyCache()
			products = rxn.GetProducts()[0]
			products.UpdatePropertyCache()

			# Assign dummy atom map numbers to unassigned reactant atoms
			current_atom_map_num = 1 + max([int(a.GetProp('molAtomMapNumber')) for a in reactants.GetAtoms() if a.HasProp('molAtomMapNumber')])
			for a in reactants.GetAtoms():
				if not a.HasProp('molAtomMapNumber'): 
					a.SetProp('molAtomMapNumber', str(current_atom_map_num))
					current_atom_map_num += 1
			reactant_smarts = Chem.MolToSmarts(reactants)

			if v: print(rxn_string)
			if v: print(template['reaction_smarts'])
			edits = summarize_reaction_outcome(reactants, products)
			if v: print(edits)
			if v: print(reactant_smarts)
						
			TRANSFORM_DB.update(
				{'_id': template['_id']},
				{
					'$set': {
						'edits': edits,
						'reactants': reactant_smarts,
						'num_edits': sum([len(z) for z in edits]),
					}
				}

			)
		except KeyboardInterrupt:
			print('User interrupted')
			break
		except Exception as e:
			print(e)

		if v: raw_input('Pause')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-v', type = str, default = 'n',
						help = 'Verbose printing; defaults to False')
	args = parser.parse_args()

	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['askcos_transforms']
	TRANSFORM_DB = db['lowe_refs_general_v3']

	v = args.v in ['Yes','yes','y','Y','1','t','T','true','True']
	main(TRANSFORM_DB, v = v)
