from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
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
		self.has_synth = False
		self.has_retro = False


		# # MEMORY LEAK DEBUG
		# from pympler.tracker import SummaryTracker
		# self.tracker = SummaryTracker()

	def load(self, collection, mincount = 4, get_retro = True, get_synth = True):
		'''
		Loads the object from a MongoDB collection containing transform
		template records.
		'''
		# Save collection source
		self.source = collection

		# Save get_retro/get_synth:
		if get_retro: self.has_retro = True
		if get_synth: self.has_synth = True

		if mincount and 'count' in collection.find_one(): 
			filter_dict = {'count': { '$gte': mincount}}
		else: 
			filter_dict = {}

		# Look for all templates in collection
		for document in collection.find(filter_dict):
			# Skip if no reaction SMARTS
			if 'reaction_smarts' not in document: continue
			reaction_smarts = str(document['reaction_smarts'])
			if not reaction_smarts: continue

			# Define dictionary
			template = {
				'name': 				document['name'] if 'name' in document else '',
				'reaction_smarts': 		reaction_smarts,
				'incompatible_groups': 	document['incompatible_groups'] if 'incompatible_groups' in document else [],
				'reference': 			document['reference'] if 'reference' in document else '',
				'references':			document['references'] if 'references' in document else [],
				'rxn_example': 			document['rxn_example'] if 'rxn_example' in document else '',
				'explicit_H': 			document['explicit_H'] if 'explicit_H' in document else False,
				'_id':	 				document['_id'] if '_id' in document else -1,
				'product_smiles':		document['product_smiles'] if 'product_smiles' in document else [],			
			}

			# Frequency/popularity score
			if 'count' in document: 
				template['count'] = document['count']
			elif 'popularity' in document:
				template['count'] = document['popularity']
			else:
				template['count'] = 1

			# Define reaction in RDKit and validate
			if get_retro:
				try:
					# Force reactants and products to be one molecule (not really, but for bookkeeping)
					reaction_smarts_retro = '(' + reaction_smarts.replace('>>', ')>>(') + ')'
					rxn = AllChem.ReactionFromSmarts(str(reaction_smarts_retro))
					#if rxn.Validate() == (0, 0):
					if rxn.Validate()[1] == 0: 
						template['rxn'] = rxn
					else:
						template['rxn'] = None
				except Exception as e:
					print('Couldnt load retro: {}: {}'.format(reaction_smarts_retro, e))
					template['rxn'] = None

			# Define forward version, too
			if get_synth:
				try:
					products, reactants = reaction_smarts.split('>>')
					reaction_smarts_forward = '(' + reactants + ')>>(' + products + ')'
					rxn_f = AllChem.ReactionFromSmarts(reaction_smarts_forward)
					#if rxn_f.Validate() == (0, 0):
					if rxn_f.Validate()[1] == 0:
						template['rxn_f'] = rxn_f
					else:
						template['rxn_f'] = None
				except Exception as e:
					print('Couldnt load forward: {}: {}'.format(reaction_smarts_forward, e))
					template['rxn_f'] = None

			# Need to have either a retro or forward reaction be valid
			if get_retro and get_synth:
				if not template['rxn'] and not template['rxn_f']: continue
			elif get_retro:
				if not template['rxn']: continue
			elif get_synth: 
				if not template['rxn_f']: continue
			else:
				raise ValueError('Cannot run Transformer.load() with get_retro = get_synth = False')

			# Add to list
			self.templates.append(template)
		self.num_templates = len(self.templates)

	def reorder(self):
		'''
		Re-orders the list of templates (self.templates) according to 
		field 'count' in descending order. This means we will apply the
		most popular templates first
		'''
		self.templates[:] = [x for x in sorted(self.templates, key = lambda z: z['count'], reverse = True)]


	def perform_retro(self, smiles):
		'''
		Performs a one-step retrosynthesis given a SMILES string of a
		target molecule by applying each transformation template
		sequentially.
		'''

		# Define mol to operate on
		mol = Chem.MolFromSmiles(smiles)
		smiles = Chem.MolToSmiles(mol, isomericSmiles = USE_STEREOCHEMISTRY) # to canonicalize

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
				continue
			if not outcomes: continue
			for j, outcome in enumerate(outcomes):
				try:
					for x in outcome:
						x.UpdatePropertyCache()
						Chem.SanitizeMol(x)
						[a.SetProp('molAtomMapNumber', a.GetProp('old_molAtomMapNumber')) \
							for (i, a) in enumerate(x.GetAtoms()) \
							if 'old_molAtomMapNumber' in a.GetPropsAsDict()]
				except Exception as e:
					print(e)
					continue
				smiles_list = []
				for x in outcome: 
					smiles_list.extend(Chem.MolToSmiles(x, isomericSmiles = USE_STEREOCHEMISTRY).split('.'))
				precursor = RetroPrecursor(
					smiles_list = sorted(smiles_list),
					template_id = template['_id'],
					num_examples = template['count'],
				)
				if '.'.join(precursor.smiles_list) == smiles: continue # no transformation
				result.add_precursor(precursor)

		return result

	def perform_forward(self, smiles, stop_if = None, progbar = False, singleonly = False):
		'''
		Performs a forward synthesis (i.e., reaction enumeration) given
		a SMILES string by applying each transformation template in 
		reverse sequentially

		stop_if - can be used for testing product matching based on 
		if the isomericSmiles matches with one of the products. It terminates
		early instead of going through all of the templates and returns True.
		'''

		# Define pseudo-molecule (single molecule) to operate on
		mol = Chem.MolFromSmiles(smiles)
		smiles = '.'.join(sorted(Chem.MolToSmiles(mol, isomericSmiles = USE_STEREOCHEMISTRY).split('.')))

		# Initialize results object
		result = ForwardResult(smiles)

		# Draw?
		if progbar:
			from tqdm import tqdm
			generator = tqdm(self.templates)
		else:
			generator = self.templates

		# Try each in turn
		for template in generator:
			# Perform
			try:
				outcomes = template['rxn_f'].RunReactants([mol])
			except Exception as e:
				#print('Forward warning: {}'.format(e))
				continue
				#print('Retro version of reaction: {}'.format(template['reaction_smarts']))
			if not outcomes: continue
			for j, outcome in enumerate(outcomes):
				try:
					for x in outcome:
						x.UpdatePropertyCache()
						Chem.SanitizeMol(x)
						[a.SetProp('molAtomMapNumber', a.GetProp('old_molAtomMapNumber')) \
							for (i, a) in enumerate(x.GetAtoms()) \
							if 'old_molAtomMapNumber' in a.GetPropsAsDict()]
				except Exception as e:
					#print(e)
					continue
				smiles_list = []
				for x in outcome: 
					smiles_list.extend(Chem.MolToSmiles(x, isomericSmiles = USE_STEREOCHEMISTRY).split('.'))
				# Reduce to largest (longest) product only?
				if singleonly: smiles_list = [max(smiles_list, key = len)]
				product = ForwardProduct(
					smiles_list = sorted(smiles_list),
					template_id = template['_id'],
					num_examples = template['count'],
				)
				if '.'.join(product.smiles_list) == smiles: continue # no transformation
				
				# Early termination?
				if stop_if:
					if stop_if in product.smiles_list: 
						print('Found true product - skipping remaining templates to apply')
						return True
				# If not early termination, we want to keep all products
				else:
					result.add_product(product)
		# Were we trying to stop early?
		if stop_if: 
			return False
		# Otherwise, return the full list of products
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
				# Just add this template_id and score
				old_product.template_ids = list(set(old_product.template_ids + 
													product.template_ids))
				old_product.num_examples += product.num_examples
				return
		# New!
		self.products.append(product)

	def return_top(self, n = 50):
		'''
		Returns the top n products as a list of dictionaries, 
		sorted by descending score
		'''
		top = []
		for (i, product) in enumerate(sorted(self.products, key = lambda x: x.num_examples, reverse = True)):
			top.append({
				'rank': i + 1,
				'smiles': '.'.join(product.smiles_list),
				'smiles_split': product.smiles_list,
				'num_examples': product.num_examples,
				'tforms': product.template_ids,
				})
			if i + 1 == n: 
				break
		return top

class ForwardProduct:
	'''
	A class to store a single forward product for reaction enumeration
	'''
	def __init__(self, smiles_list = [], template_id = -1, num_examples = 0):
		self.smiles_list = smiles_list
		self.template_ids = [template_id]
		self.num_examples = num_examples

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
				old_precursor.num_examples += precursor.num_examples
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
		for (i, precursor) in enumerate(sorted(self.precursors, \
				key = lambda x: (x.retroscore, x.num_examples), reverse = True)):
			top.append({
				'rank': i + 1,
				'smiles': '.'.join(precursor.smiles_list),
				'smiles_split': precursor.smiles_list,
				'score': precursor.retroscore,
				'num_examples': precursor.num_examples,
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
	def __init__(self, smiles_list = [], template_id = -1, num_examples = 0):
		self.retroscore = 0
		self.num_examples = num_examples
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
								- 2.00 * np.sum(np.power(ring_atoms, 1.5)) \
								- 0.00 * np.sum(np.power(chiral_centers, 2.0))