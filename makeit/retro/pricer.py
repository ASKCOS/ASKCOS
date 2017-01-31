from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import rdkit.Chem as Chem          
from rdkit.Chem import AllChem
import numpy as np
from collections import defaultdict
from tqdm import tqdm

class Pricer:
	'''
	The Pricer class is used to look up the ppg of chemicals if they
	are buyable.
	'''

	def __init__(self, by_xrn = False):

		self.chemicals_source = None
		self.buyables_source = None
		self.by_xrn = False

		self.prices = defaultdict(float) # default 0 ppg means not buyable
		self.prices_by_xrn = defaultdict(float)

	def load(self, chemicals_collection, buyables_collection):
		'''
		Loads the object from a MongoDB collection containing transform
		template records.
		'''
		# Save collection source
		self.chemicals_source = chemicals_collection
		self.buyables_source = buyables_collection


		buyable_dict = {}
		# First pull buyables source (smaller)
		for buyable_doc in tqdm(self.buyables_source.find({}, ['ppg', 'smiles', 'smiles_flat'], no_cursor_timeout = True)):
			if USE_STEREOCHEMISTRY:
				smiles = buyable_doc['smiles']
			else:
				smiles = buyable_doc['smiles_flat']
			buyable_dict[buyable_doc['_id']] = buyable_doc['ppg']
			self.prices[smiles] = buyable_doc['ppg']

		if self.by_xrn:
			# Then pull chemicals source for XRNs (larger)
			for chemical_doc in tqdm(self.chemicals_source.find({'buyable_id': {'$gt': -1}}, ['buyable_id'], no_cursor_timeout = True)):
				if 'buyable_id' not in chemical_doc: continue
				self.prices_by_xrn[chemical_doc['_id']] = buyable_dict[chemical_doc['buyable_id']]

	def lookup_smiles(self, smiles, alreadyCanonical = False):
		'''
		Looks up a price by SMILES. Tries it as-entered and then 
		re-canonicalizes it in RDKit unl ess the user specifies that
		the string is definitely already canonical.
		'''
		ppg = self.prices[smiles]
		if not alreadyCanonical:
			if not ppg:
				mol = Chem.MolFromSmiles(smiles)
				if not mol: return ppg
				smiles = Chem.MolToSmiles(mol, isomericSmiles = USE_STEREOCHEMISTRY)
				ppg = self.prices[smiles]
		return ppg

	def lookup_xrn(self, xrn):
		'''
		Looks up a price by Reaxys XRN.
		'''
		if not self.by_rxn: raise ValueError('Not initialized to look up prices by XRN!')
		return self.prices_by_xrn[xrn]