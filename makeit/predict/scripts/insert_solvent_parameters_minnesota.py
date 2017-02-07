# Import relevant packages
import numpy as np     	      	   # for simple calculations
import rdkit.Chem as Chem 
import os
import cirpy

# DATABASE
from pymongo import MongoClient    # mongodb plugin
client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaxys']
SOLVENT_DB = db['solvents_mn']

import csv
with open(os.path.join(os.path.dirname(__file__), 'Minnesota solvent parameters.csv'), 'rb') as csvfile:
	reader = csv.reader(csvfile)
	for i, row in enumerate(reader):
		if i == 0: continue # headers

		name = row[0]
		n = float(row[1])
		n25 = float(row[2]) if row[2] else float(row[1])
		a = float(row[3])
		b = float(row[4])
		g = float(row[5])
		e = float(row[6])
		phi = float(row[7])
		psi = float(row[8])

		n_norm = float(row[9])
		n25_norm = float(row[10]) if row[2] else float(row[9])
		a_norm = float(row[11])
		b_norm = float(row[12])
		g_norm = float(row[13])
		e_norm = float(row[14])
		phi_norm = float(row[15])
		psi_norm = float(row[16])
		

		smiles = None 

		# Use RDKit to canonicalize
		try:
			smiles = cirpy.resolve(name, 'smiles')
		except Exception as e:
			print(e)
			print(name)
			raw_input('Pausing...') 

		if not smiles: 
			print('## cirpy could not resolve {}'.format(name))
		mol = Chem.MolFromSmiles(smiles)
		smiles = Chem.MolToSmiles(mol)

		doc = {
			'_id'   : smiles,
			'name'  : name,
			'smiles': smiles,
			'n'     : n,
			'n25'   : n25,
			'a'     : a,
			'b'     : b,
			'g'     : g,
			'e'     : e,
			'phi'   : phi,
			'psi'   : psi,
			'n_norm'     : n_norm,
			'n25_norm'   : n25_norm,
			'a_norm'     : a_norm,
			'b_norm'     : b_norm,
			'g_norm'     : g_norm,
			'e_norm'     : e_norm,
			'phi_norm'   : phi_norm,
			'psi_norm'   : psi_norm,
		}

		try:
			SOLVENT_DB.insert_one(doc)
		except Exception as e:
			print(e)
			print(name)


print('Loaded {} solvents'.format(SOLVENT_DB.count()))