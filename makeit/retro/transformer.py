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

	def load(self, collection, mincount = 2):
		'''
		Loads the object from a MongoDB collection containing transform
		template records.
		'''
		# Save collection source
		self.source = collection

		if mincount and 'count' in collection.find_one(): 
			filter_dict = {'count': { '$gt': 2}}
		else: 
			filter_dict = {}

		# Look for all templates in collection
		for document in collection.find(filter_dict):
			# Skip if no reaction SMARTS
			if 'reaction_smarts' not in document: continue
			reaction_smarts = document['reaction_smarts']
			if not reaction_smarts: continue

			# Define dictionary
			template = {
				'name': 				document['name'] if 'name' in document else '',
				'reaction_smarts': 		reaction_smarts,
				'incompatible_groups': 	document['incompatible_groups'] if 'incompatible_groups' in document else [],
				'reference': 			document['reference'] if 'reference' in document else '',
				'rxn_example': 			document['rxn_example'] if 'rxn_example' in document else '',
				'explicit_H': 			document['explicit_H'] if 'explicit_H' in document else False,
				'_id':	 				document['_id'] if '_id' in document else -1,
				'product_smiles':		document['product_smiles'] if 'product_smiles' in document else [],
			}

			if 'count' in document: template['count'] = document['count']

			# Define reaction in RDKit and validate
			try:
				# Force products to be one molecule (not really, but for bookkeeping)
				if '.' in reaction_smarts.split('>>')[1]:
					reaction_smarts = reaction_smarts.replace('>>', '>>(') + ')'
				rxn = AllChem.ReactionFromSmarts(str(reaction_smarts))
			except:
				continue
			if rxn.Validate() != (0, 0): continue
			template['rxn'] = rxn

			# Add to list
			self.templates.append(template)
		self.num_templates = len(self.templates)

	def perform_retro(self, smiles):
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
				byproducts = [Chem.MolFromSmiles(x) for x in template['product_smiles']]
				outcomes = template['rxn'].RunReactants([mol] + byproducts)
			except Exception as e:
				print('warning: {}'.format(e))
				print(template['reaction_smarts'])
			if not outcomes: continue
			for j, outcome in enumerate(outcomes):
				try:
					[x.UpdatePropertyCache() for x in outcome]
					[Chem.SanitizeMol(x) for x in outcome]
				except Exception as e:
					print(e)
					continue
				smiles_list = []
				for x in outcome: 
					smiles_list.extend(Chem.MolToSmiles(x, isomericSmiles = True).split('.'))
				precursor = RetroPrecursor(
					smiles_list = sorted(smiles_list),
					template_id = template['_id']
				)
				if '.'.join(precursor.smiles_list) == smiles: continue # no transformation
				result.add_precursor(precursor)

		return result

	def perform_forward(self, smiles):
		'''
		Performs a forward synthesis (i.e., reaction enumeration) given
		a SMILES string by applying each transformation template in 
		reverse sequentially
		'''

		# Define pseudo-molecule (single molecule) to operate on
		mol = Chem.MolFromSmiles(smiles)
		smiles = '.'.join(sorted(Chem.MolToSmiles(mol, isomericSmiles = True).split('.')))

		# Initialize results object
		result = ForwardResult(smiles)

		# Try each in turn
		for template in self.templates:
			# Need to generate reaction
			# only retrosynthesis direction is saved as RDKit ChemicalReaction
			products, reactants = template['reaction_smarts'].split('>>')
			reaction_smarts = '(' + reactants + ')>>(' + products + ')'
			# Define reaction in RDKit and validate
			try:
				rxn = AllChem.ReactionFromSmarts(reaction_smarts)
			except Exception as e:
				print('Could not parse forward reaction: {}'.format(reaction_smarts))
				continue
			if rxn.Validate() != (0, 0): continue

			# Perform
			try:
				outcomes = rxn.RunReactants([mol])
			except Exception as e:
				print('warning: {}'.format(e))
				print(reaction_smarts)
			if not outcomes: continue
			for j, outcome in enumerate(outcomes):
				try:
					[x.UpdatePropertyCache() for x in outcome]
					[Chem.SanitizeMol(x) for x in outcome]
				except Exception as e:
					print(e)
					continue
				smiles_list = []
				for x in outcome: 
					smiles_list.extend(Chem.MolToSmiles(x, isomericSmiles = True).split('.'))
				product = ForwardProduct(
					smiles_list = sorted(smiles_list),
					template_id = template['_id']
				)
				if '.'.join(product.smiles_list) == smiles: continue # no transformation
				result.add_product(product)
		
		return result

	def lookup_id(self, template_id):
		'''
		Find the reaction smarts for this template_id
		'''
		for template in self.templates:
			if template['_id'] == template_id:
				return template

class ForwardResult:
	'''
	A class to store the results of a one-step forward synthesis.
	'''

	def __init__(self, smiles):
		self.smiles = smiles 
		self.products = []

	def add_product(self, product):
		'''
		Adds a product to the product set if it is a new product
		'''
		# Check if it is new or old
		for old_product in self.products:
			if product.smiles_list == old_product.smiles_list:
				# Just add this template_id
				old_product.template_ids = list(set(old_product.template_ids + 
													product.template_ids))
				return
		# New!
		self.products.append(product)

	def return_top(self, n = 50):
		'''
		Returns the top n products as a list of dictionaries, 
		sorted by descending score
		'''
		top = []
		np.random.shuffle(self.products)
		for (i, product) in enumerate(self.products):
			top.append({
				'rank': i + 1,
				'smiles': '.'.join(product.smiles_list),
				'smiles_split': product.smiles_list,
				'score': 0,
				'tforms': product.template_ids,
				})
			if i + 1 == n: 
				break
		return top

class ForwardProduct:
	'''
	A class to store a single forward product for reaction enumeration
	'''
	def __init__(self, smiles_list = [], template_id = -1):
		self.smiles_list = smiles_list
		self.template_ids = [template_id]

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