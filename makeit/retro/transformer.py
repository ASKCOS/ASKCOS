import rdkit.Chem as Chem          
from rdkit.Chem import AllChem
import numpy as np

class Transformer:
	'''
	The Transformer class defines an object which can be used to perform
	one-step retrosyntheses for a given molecule.
	'''

	def __init__(self):
		self.source = None
		self.templates = []

	def load(self, collection):
		'''
		Loads the object from a MongoDB collection containing transform
		template records.
		'''
		# Save collection source
		self.source = collection

		# Look for all templates in collection
		for document in collection.find():
			# Skip if no reaction SMARTS
			if 'reaction_smarts' not in document: continue
			reaction_smarts = document['reaction_smarts'] if 'reaction_smarts' in document else ''
			if not reaction_smarts: continue

			for force_intramolecular in [False, True]:
				# Force grouped products? Only if multiple products!
				if force_intramolecular: 
					if '.' in reaction_smarts.split('>>')[1]:
						reaction_smarts = reaction_smarts.replace('>>', '>>(') + ')'
					else:
						pass

				# Define dictionary
				template = {
					'name': 				document['name'] if 'name' in document else '',
					'reaction_smarts': 		reaction_smarts,
					'incompatible_groups': 	document['incompatible_groups'] if 'incompatible_groups' in document else [],
					'reference': 			document['reference'] if 'reference' in document else '',
					'rxn_example': 			document['rxn_example'] if 'rxn_example' in document else '',
					'explicit_H': 			document['explicit_H'] if 'explicit_H' in document else False,
					'_id':	 				document['_id'] if '_id' in document else -1,
				}

				# Define reaction in RDKit and validate
				try:
					rxn = AllChem.ReactionFromSmarts(str(template['reaction_smarts']))
				except:
					continue
				if rxn.Validate() != (0, 0): continue
				template['rxn'] = rxn

				# Add to list
				self.templates.append(template)
		self.num_templates = len(self.templates)

	def perform(self, smiles):
		'''
		Performs a one-step retrosynthesis given a SMILES string of a
		target molecule by applying each transformation template
		sequentially.
		'''

		# Define mol to operate on
		mol = Chem.MolFromSmiles(smiles)
		smiles = Chem.MolToSmiles(mol, isomericSmiles = True) # to canonicalize

		# Initialize results object
		result = RetroResult(smiles)

		# Try each in turn
		for template in self.templates:
			try:
				outcomes = template['rxn'].RunReactants([mol])
			except Exception as e:
				print('warning: {}'.format(e))
			if not outcomes: continue
			for j, outcome in enumerate(outcomes):
				try:
					[x.UpdatePropertyCache() for x in outcome]
					[Chem.SanitizeMol(x) for x in outcome]
				except:
					continue
				precursor = RetroPrecursor(
					smiles_list = sorted([Chem.MolToSmiles(x, isomericSmiles = True) for x in outcome]),
					template_id = template['_id']
				)
				if '.'.join(precursor.smiles_list) == smiles: continue # no transformation
				result.add_precursor(precursor)

		return result

	def lookup_id(self, template_id):
		'''
		Find the reaction smarts for this template_id
		'''
		for template in self.templates:
			if template['_id'] == template_id:
				return template['reaction_smarts']

class RetroResult:
	'''
	A class to store the results of a one-step retrosynthesis.
	'''
	def __init__(self, target_smiles):
		self.target_smiles = target_smiles
		self.precursors = []

	def add_precursor(self, precursor):
		'''
		Adds a precursor to the retrosynthesis result if it is a new
		and unique product
		'''
		# Check if the precursor set is new or old
		for old_precursor in self.precursors:
			if precursor.smiles_list == old_precursor.smiles_list:
				# Just need to add the fact that this template_id can make it
				old_precursor.template_ids = list(set(old_precursor.template_ids + 
					                                  precursor.template_ids))
				return
		# New! Need to score and add to list
		precursor.score()
		self.precursors.append(precursor)

	def return_top(self, n = 50):
		'''
		Returns the top n precursors as a list of dictionaries, 
		sorted by descending score
		'''
		top = []
		for (i, precursor) in enumerate(sorted(self.precursors, key = lambda x: x.retroscore, reverse = True)):
			top.append({
				'rank': i + 1,
				'smiles': '.'.join(precursor.smiles_list),
				'smiles_split': precursor.smiles_list,
				'score': precursor.retroscore,
				'tforms': precursor.template_ids,
				})
			if i + 1 == n: 
				break
		return top

class RetroPrecursor:
	'''
	A class to store a single set of precursor(s) for a retrosynthesis
	does NOT contain the target molecule information
	'''
	def __init__(self, smiles_list = [], template_id = -1):
		self.retroscore = 0
		self.smiles_list = smiles_list
		self.template_ids = [template_id]

	def score(self):
		'''
		Calculate the score of this step
		'''
		mols = [Chem.MolFromSmiles(x) for x in self.smiles_list]
		total_atoms = [x.GetNumHeavyAtoms() for x in mols]
		ring_atoms = [sum([a.IsInRing() for a in x.GetAtoms()])	for x in mols]
		chiral_centers = [len(Chem.FindMolChiralCenters(x)) for x in mols]
		self.retroscore = - 1.00 * np.sum(np.power(total_atoms, 1.5)) \
								- 5.00 * np.sum(np.power(ring_atoms, 1.5)) \
								- 0.00 * np.sum(np.power(chiral_centers, 2.0))