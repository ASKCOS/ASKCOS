from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import rdkit.Chem as Chem          
from rdkit.Chem import AllChem
import numpy as np

from makeit.retro.transformer import RetroResult, RetroPrecursor

class TransformerOnlyKnown:
	'''
	The TransformerOnlyKnown class defines an object which can be used to perform
	one-step retrosyntheses for a given molecule. Only reactions found in
	the literature are included
	'''

	def __init__(self):
		self.chemicals_source = None
		self.reactions_source = None

	def load(self, chemicals_collection, reactions_collection):
		'''
		Loads the object from a MongoDB collection containing chemicals/reactions
		'''
		# Save collection source
		self.chemicals_source = chemicals_collection
		self.reactions_source = reactions_collection

	def perform_retro(self, smiles):
		'''
		Performs a one-step retrosynthesis given a SMILES string of a
		target molecule by searching for reactions where it appears as
		a product
		'''

		# Add the target compound
		chem_doc = self.chemicals_source.find_one({'SMILES': smiles}, ['_id'])

		if not chem_doc: # try to canonicalize
			mol = Chem.MolFromSmiles(smiles)
			if not mol: return False
			smiles = Chem.MolToSmiles(mol)
			chem_doc = self.chemicals_source.find_one({'SMILES': smiles}, ['_id'])
			if not chem_doc: return False

		xrn = chem_doc['_id']
		print('Matched target to ID {}'.format(xrn))

		# Initialize results object
		result = RetroResult(smiles)
		
		# Look for all reactions where this is the product
		for rx_doc in self.reactions_source.find({'RX_PXRN': xrn, 'RX_SKW': {'$ne': 'half reaction'}}, ['RX_RXRN', 'RX_NVAR']):

			new_xrns = rx_doc['RX_RXRN']
			if not new_xrns: continue

			smiles_list = []
			for new_xrn in new_xrns:
				chem_doc = self.chemicals_source.find_one({'_id': new_xrn}, ['SMILES'])

				if not chem_doc: 
					print('Warning: could not find entry for ID {}'.format(new_xrn))
					continue
				if 'SMILES' not in chem_doc:
					print('Chemical ID {} does not have a SMILES string'.format(new_xrn))
					continue
				smiles_list.append(str(chem_doc['SMILES']))

			precursor = RetroPrecursor(
				smiles_list = sorted(smiles_list),
				template_id = rx_doc['_id'],
				num_examples = rx_doc['RX_NVAR'],
			)
			if '.'.join(precursor.smiles_list) == smiles: continue # no transformation
			result.add_precursor(precursor)

		return result