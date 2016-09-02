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
		c = row[2]
		e = row[3]
		s = row[4]
		a = row[5]
		b = row[6]
		v = row[7]
		mp = row[8]
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