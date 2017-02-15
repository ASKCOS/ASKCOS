import rdkit.Chem as Chem
from tqdm import tqdm
from pymongo import MongoClient

if __name__ == '__main__':
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaxys']
	CHEMICAL_DB = db['chemicals']

	for i, chem_doc in tqdm(enumerate(
			CHEMICAL_DB.find({'SMILES': {'$exists': True}, 'SMILES_FLAT': {'$exists': False}}, ['SMILES'], 
			no_cursor_timeout = True)
			)):

		# Manuually skipping when restarting
		if i < 7598678: continue 

		if 'SMILES' not in chem_doc: continue
		mol = Chem.MolFromSmiles(chem_doc['SMILES'])
		if not mol: continue
		smiles_flat = Chem.MolToSmiles(mol, isomericSmiles = False)
		smiles =      Chem.MolToSmiles(mol, isomericSmiles = True)

		to_set = {}

		if smiles != chem_doc['SMILES']:
			to_set['SMILES'] = smiles

		if smiles_flat != smiles and smiles_flat != chem_doc['SMILES']:
			to_set['SMILES_FLAT'] = smiles_flat 

		if to_set:
			# print('Chemical {}'.format(chem_doc['_id']))
			# print('Recorded SMILES: {}'.format(chem_doc['SMILES']))
			# print('Canonical SMILES: {}'.format(smiles))
			# print('Flat SMILES:     {}'.format(smiles_flat))
			# print('to_set: {}'.format(to_set))
			# raw_input('Pause...')
			CHEMICAL_DB.update_one(
				{'_id': chem_doc['_id']},
				{'$set': to_set},
			)

	print('{} chemicals now have a SMILES_FLAT != SMILES'.format(
		CHEMICAL_DB.count({'SMILES_FLAT': {'$exists': True}})
	))