# Import relevant packages
import numpy as np     	      	   # for simple calculations
import rdkit.Chem as Chem 
import os

# DATABASE
from pymongo import MongoClient    # mongodb plugin
client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaxys']
SOLVENT_DB = db['solvents']

import csv
with open(os.path.join(os.path.dirname(__file__), 'Abraham Solvent Parameters to load.csv'), 'rb') as csvfile:
	reader = csv.reader(csvfile)
	for i, row in enumerate(reader):
		if i == 0: continue # headers

		name = row[0]
		smiles = row[1]
		c = float(row[2])
		e = float(row[3])
		s = float(row[4])
		a = float(row[5])
		b = float(row[6])
		v = float(row[7])
		try:
			mp = float(row[8])
		except Exception:
			mp = row[8]
		try:
			bp = float(row[9])
		except Exception:
			bp = row[9]

		# Use RDKit to canonicalize
		mol = Chem.MolFromSmiles(smiles)
		smiles = Chem.MolToSmiles(mol)

		doc = {
			'_id': smiles,
			'name': name,
			'smiles': smiles,
			'c': c,
			'e': e,
			's': s,
			'a': a,
			'b': b,
			'v': v,
			'mp': mp,
			'bp': bp,
		}

		try:
			SOLVENT_DB.insert_one(doc)
		except Exception as e:
			print(e)

doc = {
	'_id': 'default',
	'c': 0.2216,
	'e': 0.2677,
	's': -0.4,
	'a': -1.135,
	'b': -3.799,
	'v': 3.6212,
}
try:
	SOLVENT_DB.insert_one(doc)
except Exception as e:
	print(e)

print('Loaded {} solvents (including one default)'.format(SOLVENT_DB.count()))