from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import sys
import os
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import re
import itertools


if __name__ == '__main__':
	out_folder = os.path.join(os.path.dirname(__file__), 'output')

	# DATABASE
	from pymongo import MongoClient    # mongodb plugin
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaxys']
	TRANSFORM_DB = db['transforms_forward_v1']
	INSTANCE_DB = db['instances']
	CHEMICAL_DB = db['chemicals']
	REDUCTANT_DB = db['reductants']
	done_references = [doc['references'] for doc in REDUCTANT_DB.find({}, ['references'])]
	done_references = list(itertools.chain.from_iterable(done_references))

	Hs_change = []
	for template in TRANSFORM_DB.find(no_cursor_timeout=True):
		# print(str(template['reaction_smarts']))

		### Use RDKit to count hydrogens

		# rxn = AllChem.ReactionFromSmarts(str(template['reaction_smarts']))

		# [mol.UpdatePropertyCache() for mol in rxn.GetReactants()]
		# [mol.UpdatePropertyCache() for mol in rxn.GetProducts()]
		
		# # Get number of hydrogens on non-leaving groups
		# Hs_before = sum([sum([a.GetTotalNumHs() for a in mol.GetAtoms() if a.HasProp('molAtomMapNumber')]) for mol in rxn.GetReactants()])
		# Hs_after  = sum([sum([a.GetTotalNumHs() for a in mol.GetAtoms() if a.HasProp('molAtomMapNumber')]) for mol in rxn.GetProducts()])

		### Count number of hydrogens on atom mapped atoms
		# *except* don't count OH groups, because hydrolysis is boring
		def match_to_num_Hs(match):
			if match == 'H': return 1
			if 'O' in match: return 0
			return int(match[1:])

		Hs_before = sum([
			match_to_num_Hs(match)
			for match in re.findall('(O?H[0-9]*).*?\:', str(template['reaction_smarts']).split('>')[0], re.DOTALL)
		])
		Hs_after = sum([
			match_to_num_Hs(match)
			for match in re.findall('(O?H[0-9]*).*?\:', str(template['reaction_smarts']).split('>')[2], re.DOTALL)
		])
		Hs_change.append(Hs_after - Hs_before)

		if (Hs_after - Hs_before) >= 1:
			print('### REDUCTION ### {}->{} ### - {}'.format(Hs_before, Hs_after, template['reaction_smarts']))
			print('  processing {} references...'.format(len(template['references'])))
			for rxd_id in template['references']:
				if rxd_id in done_references: continue # save time if already done
				rxd = INSTANCE_DB.find_one({'_id': rxd_id})
				if not rxd: continue

				for xrn in rxd['RXD_RGTXRN']:
					doc = REDUCTANT_DB.find_one({'_id': xrn})
					if not doc:
						new_doc = CHEMICAL_DB.find_one({'_id': xrn}, ['_id', 'SMILES'])
						if not new_doc: continue
						new_doc['references'] = [rxd_id]
						new_doc['count'] = 1
						REDUCTANT_DB.insert(new_doc)
					elif rxd_id not in doc['references']:
						REDUCTANT_DB.update_one(
							{'_id': xrn},
							{
								'$inc': {
									'count': 1,
								},
								'$addToSet': {
									'references': rxd_id,
								},
							}
						)