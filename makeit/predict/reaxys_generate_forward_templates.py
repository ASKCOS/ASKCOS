'''
This script (generate_reaction_templates) uses a MongoDB collection 
with atom-mapped reaction SMILES strings and parses them into a new 
collection containing the transforms.
'''

from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import argparse
from numpy.random import shuffle # for random selection
import rdkit.Chem as Chem          # molecule building
from rdkit.Chem import AllChem
from collections import defaultdict
import rdkit.Chem.Draw as Draw
from rdkit import RDLogger
import datetime # for info files
import json # for dumping
import sys  # for commanad line
import os   # for file paths
import re 
import itertools
from makeit.retro.draw import *

# DATABASE
from pymongo import MongoClient    # mongodb plugin
client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaxys']
TRANSFORM_DB = db['transforms_forward_v1']
REACTION_DB = db['reactions']
INSTANCE_DB = db['instances']
CHEMICAL_DB = db['chemicals']
reaction_smiles_field = 'RXN_SMILES'

def mols_from_smiles_list(all_smiles):
	'''Given a list of smiles strings, this function creates rdkit
	molecules'''
	mols = []
	for smiles in all_smiles:
		if not smiles: continue
		mols.append(Chem.MolFromSmiles(smiles))
	return mols

def bond_to_label(bond):
	'''This function takes an RDKit bond and creates a label describing
	the most important attributes'''
	# atoms = sorted([atom_to_label(bond.GetBeginAtom().GetIdx()), \
	# 			    atom_to_label(bond.GetEndAtom().GetIdx())])
	a1_label = str(bond.GetBeginAtom().GetAtomicNum())
	a2_label = str(bond.GetEndAtom().GetAtomicNum())
	if bond.GetBeginAtom().HasProp('molAtomMapNumber'):
		a1_label += bond.GetBeginAtom().GetProp('molAtomMapNumber')
	if bond.GetEndAtom().HasProp('molAtomMapNumber'):
		a2_label += bond.GetEndAtom().GetProp('molAtomMapNumber')
	atoms = sorted([a1_label, a2_label])

	return '{}{}{}'.format(atoms[0], bond.GetSmarts(), atoms[1])

def get_tagged_atoms_from_mols(mols):
	'''Takes a list of RDKit molecules and returns total list of
	atoms and their tags'''
	atoms = []
	atom_tags = []
	for mol in mols:
		new_atoms, new_atom_tags = get_tagged_atoms_from_mol(mol)
		atoms += new_atoms 
		atom_tags += new_atom_tags
	return atoms, atom_tags

def get_tagged_atoms_from_mol(mol):
	'''Takes an RDKit molecule and returns list of tagged atoms and their
	corresponding numbers'''
	atoms = []
	atom_tags = []
	for atom in mol.GetAtoms():
		if ':' in atom.GetSmarts():
			atoms.append(atom)
			atom_tags.append(atom.GetSmarts().split(':')[1][:-1])
	return atoms, atom_tags

def atoms_are_different(atom1, atom2, level = 1):
	'''Compares two RDKit atoms based on basic properties'''

	if atom1.GetSmarts() != atom2.GetSmarts(): return True # should be very general
	if atom1.GetAtomicNum() != atom2.GetAtomicNum(): return True # must be true for atom mapping
	if atom1.GetTotalNumHs() != atom2.GetTotalNumHs(): return True
	if atom1.GetFormalCharge() != atom2.GetFormalCharge(): return True
	if atom1.GetDegree() != atom2.GetDegree(): return True
	#if atom1.IsInRing() != atom2.IsInRing(): return True
	if atom1.GetNumRadicalElectrons() != atom2.GetNumRadicalElectrons(): return True
	# TODO: add # pi electrons like ICSynth?

	# Check bonds and nearest neighbor identity
	if level >= 1:
		bonds1 = sorted([bond_to_label(bond) for bond in atom1.GetBonds()]) 
		bonds2 = sorted([bond_to_label(bond) for bond in atom2.GetBonds()]) 
		if bonds1 != bonds2: return True

		# # Check neighbors too (already taken care of with previous lines)
		# neighbors1 = sorted([atom.GetAtomicNum() for atom in atom1.GetNeighbors()])
		# neighbors2 = sorted([atom.GetAtomicNum() for atom in atom2.GetNeighbors()])
		# if neighbors1 != neighbors2: return True

	return False

def get_changed_atoms(reactants, products):
	'''Looks at mapped atoms in a reaction and determines which ones changed'''

	err = 0
	prod_atoms, prod_atom_tags = get_tagged_atoms_from_mols(products)

	if v: print('Products contain {} tagged atoms'.format(len(prod_atoms)))
	if v: print('Products contain {} unique atom numbers'.format(len(set(prod_atom_tags))))

	reac_atoms, reac_atom_tags = get_tagged_atoms_from_mols(reactants)
	if len(set(prod_atom_tags)) != len(set(reac_atom_tags)):
		if v: print('warning: different atom tags appear in reactants and products')
		#err = 1 # okay for Reaxys, since Reaxys creates mass
	if len(prod_atoms) != len(reac_atoms):
		if v: print('warning: total number of tagged atoms differ, stoichometry != 1?')
		#err = 1

	# Find differences 
	changed_atoms = []
	changed_atom_tags = []
	#print(reac_atom_tags)
	#print(prod_atom_tags)

	# Product atoms that are different from reactant atom equivalent
	for i, prod_tag in enumerate(prod_atom_tags):

		for j, reac_tag in enumerate(reac_atom_tags):
			if reac_tag != prod_tag: continue
			if reac_tag not in changed_atom_tags: # don't bother comparing if we know this atom changes
				# If atom changed, add
				if atoms_are_different(prod_atoms[i], reac_atoms[j]):
					changed_atoms.append(reac_atoms[j])
					changed_atom_tags.append(reac_tag)
					break
				# If reac_tag appears multiple times, add (need for stoichometry > 1)
				if prod_atom_tags.count(reac_tag) > 1:
					changed_atoms.append(reac_atoms[j])
					changed_atom_tags.append(reac_tag)
					break

	# Reactant atoms that do not appear in product (tagged leaving groups)
	for j, reac_tag in enumerate(reac_atom_tags):
		if reac_tag not in changed_atom_tags:
			if reac_tag not in prod_atom_tags:
				changed_atoms.append(reac_atoms[j])
				changed_atom_tags.append(reac_tag)

	if v: 
		print('{} tagged atoms in reactants change 1-atom properties'.format(len(changed_atom_tags)))
		for smarts in [atom.GetSmarts() for atom in changed_atoms]:
			print('  {}'.format(smarts))

	return changed_atoms, changed_atom_tags, err

def expand_atoms_to_use(mol, atoms_to_use, groups = [], symbol_replacements = []):
	'''Given an RDKit molecule and a list of AtomIdX which should be included
	in the reaction, this function expands the list of AtomIdXs to include one 
	nearest neighbor with special consideration of (a) unimportant neighbors and
	(b) important functional groupings'''

	# Copy
	new_atoms_to_use = atoms_to_use[:]

	# Look for all atoms in the current list of atoms to use
	for atom in mol.GetAtoms():
		if atom.GetIdx() not in atoms_to_use: continue
		# Look for all nearest neighbors of the currently-included atoms
		for neighbor in atom.GetNeighbors():
			# Evaluate nearest neighbor atom to determine what should be included
			new_atoms_to_use, symbol_replacements = \
					expand_atoms_to_use_atom(mol, new_atoms_to_use, neighbor.GetIdx(), 
						groups = groups, symbol_replacements = symbol_replacements)
			
	return new_atoms_to_use, symbol_replacements

def expand_atoms_to_use_atom(mol, atoms_to_use, atom_idx, groups = [], symbol_replacements = []):
	'''Given an RDKit molecule and a list of AtomIdx which should be included
	in the reaction, this function extends the list of atoms_to_use by considering 
	a candidate atom extension, atom_idx'''

	# Skip current candidate atom if it is already included
	if atom_idx in atoms_to_use:
		return atoms_to_use, symbol_replacements

	# See if this atom belongs to any groups
	found_in_group = False
	for group in groups:
		if int(atom_idx) in group: # int correction
			if v: print('added group centered at {}'.format(atom_idx))
			# Add the whole list, redundancies don't matter
			atoms_to_use.extend(list(group)) 
			found_in_group = True
	if found_in_group:	return atoms_to_use, symbol_replacements
		
	# How do we add an atom that wasn't in an identified important functional group?
	# Develop special SMARTS symbol

	# Include this atom
	atoms_to_use.append(atom_idx)

	# Look for replacements
	symbol_replacements.append((atom_idx, convert_atom_to_wildcard(mol.GetAtomWithIdx(atom_idx))))

	return atoms_to_use, symbol_replacements

def convert_atom_to_wildcard(atom):
	'''This function takes an RDKit atom and turns it into a wildcard 
	using hard-coded generalization rules. This function should be used
	when candidate atoms are used to extend the reaction core for higher
	generalizability'''

	if SUPER_GENERAL_TEMPLATES:
		if ':' in atom.GetSmarts():
			# Fully generalize atom-mappped neighbors because they are aren't a leaving group
			label = re.search('\:[0-9]+\]', atom.GetSmarts())
			return '[*{}'.format(label.group())

	if not SUPER_GENERAL_TEMPLATES:

		# Is this a terminal atom? We can tell if the degree is one
		if atom.GetDegree() == 1:
			return atom.GetSmarts()

	# Initialize
	symbol = '['

	# Add atom primitive (don't use COMPLETE wildcards)
	if atom.GetAtomicNum() != 6:
		symbol += '#{};'.format(atom.GetAtomicNum())
	elif atom.GetIsAromatic():
		symbol += 'c;'
	else:
		symbol += 'C;'

	if not SUPER_GENERAL_TEMPLATES:
		# Charge is important
		if atom.GetFormalCharge() != 0:
			charges = re.search('([-+]+[1-9]?)', atom.GetSmarts())
			if charges: symbol += charges.group() + ';'

	# Strip extra semicolon
	if symbol[-1] == ';': symbol = symbol[:-1]

	# Close with label or with bracket
	label = re.search('\:[0-9]+\]', atom.GetSmarts())
	if label: 
		symbol += label.group()
	else:
		symbol += ']'

	if v: 
		if symbol != atom.GetSmarts():
			print('Improved generality of atom SMARTS {} -> {}'.format(atom.GetSmarts(), symbol))

	return symbol

def get_fragments_for_changed_atoms(mols, changed_atom_tags, radius = 0, 
	category = 'reactants', expansion = [], product_new_atoms = []):
	'''Given a list of RDKit mols and a list of changed atom tags, this function
	computes the SMILES string of molecular fragments using MolFragmentToSmiles 
	for all changed fragments.

	expansion: atoms added during reactant expansion that should be included and
	           generalized in product fragment
	product_new_atoms: atoms contributed by reagents, which are handled specially
	           and should not be generalized in the product fragment'''

	fragments = ''
	for mol in mols:
		# Initialize list of replacement symbols (updated during expansion)
		symbol_replacements = []

		# Are we looking for groups? (reactants only)
		if category == 'reactants':
			groups = get_special_groups(mol)
		else:
			groups = []

		# Build list of atoms to use
		atoms_to_use = []
		for atom in mol.GetAtoms():
			# Check self (only tagged atoms)
			if ':' in atom.GetSmarts():
				if atom.GetSmarts().split(':')[1][:-1] in changed_atom_tags + product_new_atoms:
					atoms_to_use.append(atom.GetIdx())
					symbol = atom.GetSmarts()
					# CUSTOM SYMBOL CHANGES
					if atom.GetTotalNumHs() == 0:
						# Be explicit when there are no hydrogens
						if ':' in symbol: # stick H0 before label
							symbol = symbol.replace(':', ';H0:')
						else: # stick before end
							symbol = symbol.replace(']', ';H0]')
						
						# print('Being explicit about H0!!!!')
					if atom.GetFormalCharge() == 0:
						# Also be explicit when there is no charge
						if ':' in symbol: 
							symbol = symbol.replace(':', ';+0:')
						else:
							symbol = symbol.replace(']', ';+0]')
					if USE_STEREOCHEMISTRY:
						if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
							# Be explicit when there is a tetrahedral chiral tag
							if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
								tag = '@'
							elif atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
								tag = '@@'
							if ':' in symbol:
								symbol = symbol.replace(':', ';{}:'.format(tag))
							else:
								symbol = symbol.replace(']', ';{}]'.format(tag))
					if symbol != atom.GetSmarts():
						symbol_replacements.append((atom.GetIdx(), symbol))
					continue

		# Check neighbors (any atom)
		for k in range(radius):
			atoms_to_use, symbol_replacements = expand_atoms_to_use(mol, atoms_to_use, 
				groups = groups, symbol_replacements = symbol_replacements)

		if category == 'products':
			# Add extra labels to include (for products only)
			if expansion:
				for atom in mol.GetAtoms():
					if ':' not in atom.GetSmarts(): continue
					label = atom.GetSmarts().split(':')[1][:-1]
					if label in expansion and label not in changed_atom_tags:
						atoms_to_use.append(atom.GetIdx())
						# Make the expansion a wildcard
						symbol_replacements.append((atom.GetIdx(), convert_atom_to_wildcard(atom)))	
						if v: print('expanded label {} to wildcard in products'.format(label))
			
			# Make sure unmapped atoms are included (from products)
			for atom in mol.GetAtoms():
				if not atom.HasProp('molAtomMapNumber'): 
					atoms_to_use.append(atom.GetIdx())

		# Define new symbols to replace terminal species with wildcards
		# (don't want to restrict templates too strictly)
		symbols = [atom.GetSmarts() for atom in mol.GetAtoms()]
		for (i, symbol) in symbol_replacements:
			symbols[i] = symbol

		if not atoms_to_use: continue
		# if v:
		# 	print('~~ replacement for this ' + category[:-1])
		# 	print('{} -> {}'.format([mol.GetAtomWithIdx(x).GetSmarts() for (x, s) in symbol_replacements], 
		# 		                    [s for (x, s) in symbol_replacements]))
		# Remove molAtomMapNumber before canonicalization
		[x.ClearProp('molAtomMapNumber') for x in mol.GetAtoms()]
		
		fragments += '(' + AllChem.MolFragmentToSmiles(mol, atoms_to_use, 
		atomSymbols = symbols, allHsExplicit = True, 
		isomericSmiles = USE_STEREOCHEMISTRY, allBondsExplicit = True) + ').'
		
	return fragments[:-1]

def expand_changed_atom_tags(changed_atom_tags, reactant_fragments):
	'''Given a list of changed atom tags (numbers as strings) and a string consisting
	of the reactant_fragments to include in the reaction transform, this function 
	adds any tagged atoms found in the reactant side of the template to the 
	changed_atom_tags list so that those tagged atoms are included in the products'''

	expansion = []
	atom_tags_in_reactant_fragments = re.findall('\:([[0-9]+)\]', reactant_fragments)
	for atom_tag in atom_tags_in_reactant_fragments:
		if atom_tag not in changed_atom_tags:
			expansion.append(atom_tag)
	if v: print('after building reactant fragments, additional labels included: {}'.format(expansion))
	return expansion

def get_special_groups(mol):
	'''Given an RDKit molecule, this function returns a list of tuples, where
	each tuple contains the AtomIdx's for a special group of atoms which should 
	be included in a fragment all together. This should only be done for the 
	reactants, otherwise the products might end up with mapping mismatches'''

	if SUPER_GENERAL_TEMPLATES: return []

	# Define templates, based on Functional_Group_Hierarchy.txt from Greg Laandrum
	group_templates = [ 
		'C(=O)Cl', # acid chloride
		'C(=O)[O;H,-]', # carboxylic acid
		'[$(S-!@[#6])](=O)(=O)(Cl)', # sulfonyl chloride
		'[$(B-!@[#6])](O)(O)', # boronic acid
		'[$(N-!@[#6])](=!@C=!@O)', # isocyanate
		'[N;H0;$(N-[#6]);D2]=[N;D2]=[N;D1]', # azide
		'O=C1N(Br)C(=O)CC1', # NBS brominating agent
		'C=O', # carbonyl
		'ClS(Cl)=O', # thionyl chloride
		'[Mg][Br,Cl]', # grinard (non-disassociated)
		'[#6]S(=O)(=O)[O]', # RSO3 leaving group
		'[O]S(=O)(=O)[O]', # SO4 group
		'[N-]=[N+]=[C]', # diazo-alkyl
		]

	# Build list
	groups = []
	for template in group_templates:
		matches = mol.GetSubstructMatches(Chem.MolFromSmarts(template))
		groups.extend(list(matches))
	return groups

def mol_list_to_inchi(mols):
	'''List of RDKit molecules to InChI string separated by ++'''
	inchis = [Chem.MolToInchi(mol) for mol in mols]
	return ' ++ '.join(sorted(inchis))

def mol_list_from_inchi(inchis):
	'''InChI string separated by ++ to list of RDKit molecules'''
	return [Chem.MolFromInchi(inchi.strip()) for inchi in inchis.split('++')]

def canonicalize_template(template):
	'''This function takes one-half of a template SMARTS string 
	(i.e., reactants or products) and re-orders them based on
	an equivalent string without atom mapping.'''

	# Strip labels to get sort orders
	template_nolabels = re.sub('\:[0-9]+\]', ']', template)

	# Split into separate molecules *WITHOUT wrapper parentheses*
	template_nolabels_mols = template_nolabels[1:-1].split(').(')
	template_mols          = template[1:-1].split(').(')

	# Split into fragments within those molecules
	for i in range(len(template_mols)):
		nolabel_mol_frags = template_nolabels_mols[i].split('.')
		mol_frags         = template_mols[i].split('.')

		# Get sort order within molecule, defined WITHOUT labels
		sortorder = [j[0] for j in sorted(enumerate(nolabel_mol_frags), key = lambda x:x[1])]

		# Apply sorting and merge list back into overall mol fragment
		template_nolabels_mols[i] = '.'.join([nolabel_mol_frags[j] for j in sortorder])
		template_mols[i]          = '.'.join([mol_frags[j] for j in sortorder])

	# Get sort order between molecules, defined WITHOUT labels
	sortorder = [j[0] for j in sorted(enumerate(template_nolabels_mols), key = lambda x:x[1])]

	# Apply sorting and merge list back into overall transform
	template = '(' + ').('.join([template_mols[i] for i in sortorder]) + ')'

	return template

def reassign_atom_mapping(transform):
	'''This function takes an atom-mapped reaction and reassigns 
	the atom-mapping labels (numbers) from left to right, once 
	that transform has been canonicalized.'''

	all_labels = re.findall('\:([0-9]+)\]', transform)

	# Define list of replacements which matches all_labels *IN ORDER*
	replacements = []
	replacement_dict = {}
	counter = 1
	for label in all_labels: # keep in order! this is important
		if label not in replacement_dict:
			replacement_dict[label] = str(counter)
			counter += 1
		replacements.append(replacement_dict[label])

	# Perform replacements in order
	transform_newmaps = re.sub('\:[0-9]+\]', 
		lambda match: (':' + replacements.pop(0) + ']'),
		transform)

	return transform_newmaps

def canonicalize_transform(transform):
	'''This function takes an atom-mapped SMARTS transform and
	converts it to a canonical form by, if nececssary, rearranging
	the order of reactant and product templates and reassigning
	atom maps.'''

	transform_reordered = '>>'.join([canonicalize_template(x) for x in transform.split('>>')])
	return reassign_atom_mapping(transform_reordered)

def main(N = 15):
	'''Read reactions'''
	global v

	# Define scoring variables
	total_attempted = 0 # total reactions simulated (excludes skipped)
	total_correct = 0 # actual products predicted
	total_precise = 0 # ONLY actual products predicted
	total_unmapped = 0
	total_partialmapped = 0
	total_nonreaction = 0

	N = min([N, REACTION_DB.count()])

	try: # to allow breaking
		# Look for entries
		i = -1
		for example_doc in REACTION_DB.find(no_cursor_timeout=True):
			i += 1

			# Temporary (skipping already done)
			#if i < 1512985: continue
	
			# Are we done?
			if i == N:
				i -= 1
				break

			if v: 
				print('##################################')
				print('###        RXN {}'.format(i))
				print('##################################')

			if "nonmapped reaction" in example_doc['RX_SKW']:
				print('Unmapped reaction!')
				print('ID: {}'.format(example_doc['_id']))
				total_unmapped += 1
				continue

			try:
				# Unpack
				reaction_smiles = str(example_doc[reaction_smiles_field])
				reactants, agents, products = [mols_from_smiles_list(x) for x in 
											[mols.split('.') for mols in reaction_smiles.split('>')]]
				[Chem.SanitizeMol(mol) for mol in reactants + agents + products]
			except Exception as e:
				# can't sanitize -> skip
				print(e)
				print('Could not load SMILES or sanitize')
				print('ID: {}'.format(example_doc['_id']))
				continue
			
			try:
				###
				### Check product atom mapping to see if reagent contributes
				###

				are_unmapped_product_atoms = False
				extra_reactant_fragment = ''
				for product in products:
					if sum([a.HasProp('molAtomMapNumber') for a in product.GetAtoms()]) < len(product.GetAtoms()):
						print('Not all product atoms have atom mapping')
						print('ID: {}'.format(example_doc['_id']))
						if v: print('REACTION: {}'.format(example_doc['RXN_SMILES']))
						are_unmapped_product_atoms = True
				if are_unmapped_product_atoms and not ALLOW_CREATE_MASS: # continue
					total_partialmapped += 1
					continue
				if are_unmapped_product_atoms and ALLOW_CREATE_MASS: # add fragment to template

					total_partialmapped += 1
					starting_map = 1000
					for product in products:
						# Get unmapped atoms
						unmapped_ids = [
							a.GetIdx() for a in product.GetAtoms() if not a.HasProp('molAtomMapNumber')
						]
						# Add new atom maps
						[a.SetProp('molAtomMapNumber', str(starting_map+z)) \
							for (z, a) in enumerate(product.GetAtoms()) if not a.HasProp('molAtomMapNumber')]
						starting_map += len(unmapped_ids)
						# Define new atom symbols for fragment with atom maps, generalizing fully
						atom_symbols = ['[{}:{}]'.format(a.GetSymbol(), a.GetProp('molAtomMapNumber')) 
							for (z, a) in enumerate(product.GetAtoms())]
						# And bond symbols...
						bond_symbols = ['~' for b in product.GetBonds()]
						extra_reactant_fragment += \
							AllChem.MolFragmentToSmiles(product, unmapped_ids, 
							allHsExplicit = False, isomericSmiles = USE_STEREOCHEMISTRY, 
							atomSymbols = atom_symbols, bondSymbols = bond_symbols) + '.'
					if extra_reactant_fragment:
						extra_reactant_fragment = '(' + extra_reactant_fragment[:-1] + ')'
						reaction_smiles = extra_reactant_fragment + reaction_smiles
						if v: print('    extra reactant fragment: {}'.format(extra_reactant_fragment))

					# Look for duplicated reactant fragments that might appear multiple times
					# as a result of stoichometry
					fragments = extra_reactant_fragment[1:-1].split('.')
					fragments_labels = [re.findall('\:([[0-9]+)\]', fragment) for fragment in fragments]
					fragments_nolabels = [re.sub('\:[0-9]+\]', ']', fragment) for fragment in fragments]
					replacements = []
					for i, fragments_nolabel in enumerate(fragments_nolabels):
						for j, fragments_nolabel2 in enumerate(fragments_nolabels):
							if j <= i: continue
							if fragments_nolabel == fragments_nolabel2:
								# print('Duplicate fragments! {} and {}'.format(
								# 		fragments[i], fragments[j]
								# 	))
								for k in range(len(fragments_labels[j])):
									replacements.append((fragments_labels[j][k], fragments_labels[i][k]))
					for replacement in replacements[::-1]:
						print(replacement)
						# Perform replacements **ON REACTANT FRAGMENT** from last to first
						for i, fragment in enumerate(fragments):
							# print('Initial fragment: {}'.format(fragments[i]))
							fragments[i] = fragment.replace(':{}]'.format(replacement[0]), ':{}]'.format(replacement[1]))
							# print('After replacement: {}'.format(fragments[i]))
						# Perform replacements **ON PRODUCT MAPPING**
						for product in products:
							for atom in product.GetAtoms():
								if atom.HasProp('molAtomMapNumber'):
									if str(atom.GetProp('molAtomMapNumber')) == str(replacement[0]):
										atom.SetProp('molAtomMapNumber', str(replacement[1]))
										if v: print('Replaced atom map in products {} with {}'.format(replacement[0], replacement[1]))
					# Recombine fragments, using set to eliminate identical ones
					extra_reactant_fragment = '(' + '.'.join(set(fragments)) + ')'
					if v: print('     after consolidating duplicate fragments: {}'.format(extra_reactant_fragment))
					fragmatch = Chem.MolFromSmarts(extra_reactant_fragment[1:-1]) # no parentheses

				###
				### Do RX-level processing
				###  

				# Get unique InChi for products we are looking for
				# remove cases where product shows up
				product_inchis_split = mol_list_to_inchi(products).split(' ++ ')
				for reactant in reactants:
					reactant_inchi = mol_list_to_inchi([reactant])
					if reactant_inchi in product_inchis_split:
						product_inchis_split.remove(reactant_inchi)
				product_inchis = ' ++ '.join(product_inchis_split)

				if v: print(reaction_smiles)
				if None in reactants + products:
					print('Could not parse all molecules in reaction, skipping')
					print('ID: {}'.format(example_doc['_id']))
					continue
				if not product_inchis: 
					print('Product molecules no different than reactant molecules')
					print('ID: {}'.format(example_doc['_id']))
					total_nonreaction += 1
					continue

				# Calculate changed atoms
				changed_atoms, changed_atom_tags, err = get_changed_atoms(reactants, products)
				if err: 
					print('Could not get changed atoms')
					print('ID: {}'.format(example_doc['_id']))
					continue
				if not changed_atom_tags:
					print('No atoms changed?')
					print('ID: {}'.format(example_doc['_id']))
					continue

				# Get fragments for reactants
				reactant_fragments = get_fragments_for_changed_atoms(reactants, changed_atom_tags, 
					radius = 1, expansion = [], category = 'reactants')
				# Get fragments for products 
				# (WITHOUT matching groups but WITH the addition of reactant fragments)
				# "product_new_atoms" has the previously unmapped product atoms
				product_fragments  = get_fragments_for_changed_atoms(products, changed_atom_tags, 
					radius = 0, expansion = expand_changed_atom_tags(changed_atom_tags, reactant_fragments),
					category = 'products', 
					product_new_atoms = re.findall('\:([[0-9]+)\]', extra_reactant_fragment))

				###
				### Now look for specific RXDs.
				###
				rxd_id_list = ['{}-{}'.format(example_doc['_id'], j) for j in range(1, int(example_doc['RX_NVAR']) + 1)]

				# Look up instances of this reaction to see the reagent(s)
				for rxd_id in rxd_id_list:

					# Initialize 'reagent' to we can see if we need to add one to the reactant mixture
					reagent = None

					# Only need to do this if there are reagent-contributed atoms
					extra_reactant_fragment_rxd = extra_reactant_fragment
					if extra_reactant_fragment:
						rxd = INSTANCE_DB.find_one({'_id': rxd_id})
						print(' Instance {}'.format(rxd['_id']))
						reagents = []
						num_matches = []
						for xrn in rxd['RXD_RGTXRN']:
							smi = CHEMICAL_DB.find_one({'_id': xrn})['SMILES']
							if v: print('    reagent: {}'.format(smi))
							for smi_split in smi.split('.'):
								m = Chem.MolFromSmiles(str(smi_split))
								if m: 
									reagents.append(m)
									num_match = len(m.GetSubstructMatches(fragmatch))
									num_matches.append(num_match)
								else:
									if v: print('        (could not parse {})'.format(smi_split))
						if sum(num_matches) == 0: 
							print('None of the reagents have {}'.format(extra_reactant_fragment[1:-1]))
							print('...skipping this RXD')
							continue
						for i, reagent in enumerate(reagents):
							if num_matches[i] == max(num_matches):
								if v: print('    most likely source of {} -> {}'.format(extra_reactant_fragment[1:-1], Chem.MolToSmiles(reagent)))
								# Get a match
								reagent_ids = reagent.GetSubstructMatch(fragmatch)
								# Assign atom mapping to reagent, doesn't matter much
								for j, label in enumerate(re.findall('\:([[0-9]+)\]', extra_reactant_fragment)):
									reagent.GetAtomWithIdx(reagent_ids[j]).SetProp('molAtomMapNumber', str(label))
								extra_reactant_fragment_rxd = Chem.MolToSmiles(reagent, allHsExplicit = True, 
									isomericSmiles = USE_STEREOCHEMISTRY)
								print('For this RXD, extra fragment = {}'.format(extra_reactant_fragment_rxd))
								break

					# Report transform
					if extra_reactant_fragment_rxd:
						rxn_string = '{}.{}>>{}'.format(reactant_fragments, extra_reactant_fragment_rxd, product_fragments)
					else:
						rxn_string = '{}>>{}'.format(reactant_fragments, product_fragments)
					if v: print('\nOverall fragment transform: {}'.format(rxn_string))

					# Load into RDKit
					rxn = AllChem.ReactionFromSmarts(rxn_string)
					if rxn.Validate()[1] != 0: 
						print('Could not validate reaction successfully')
						print('ID: {}'.format(example_doc['_id']))
						print('rxn_string: {}'.format(rxn_string))
						print('original: {}'.format(example_doc['RXN_SMILES']))
						if v: raw_input('Pausing...')
						continue
					# print('here')

					# # Analyze
					# if v: 
					# 	print('Original number of reactants: {}'.format(len(reactants)))
					# 	print('Number of reactants in transform: {}'.format(rxn.GetNumReactantTemplates()))

					# Run reaction
					correct = False
					if not DO_NOT_TEST:

						# Add reacting reagent?
						if reagent:
							reactants.append(reagent)
							if v: print('Added {} to the reactant mixture to test reaction'.format(Chem.MolToSmiles(reagent)))

						# Try all combinations of reactants that fit template
						combinations = itertools.combinations(reactants, rxn.GetNumReactantTemplates())
						unique_product_sets = []
						for combination in combinations:
							outcomes = rxn.RunReactants(list(combination))
							if not outcomes: continue
							#if v: print('\nFor reactants {}'.format([Chem.MolToSmiles(mol) for mol in combination]))
							for j, outcome in enumerate(outcomes):
								#if v: print('- outcome {}/{}'.format(j + 1, len(outcomes)))
								for k, product in enumerate(outcome):
									if v: print('~~prod: {}'.format(Chem.MolToSmiles(product)))
									try:
										Chem.SanitizeMol(product)
										product.UpdatePropertyCache()
									except Exception as e:
										print('warning: {}'.format(e))
									if v: 
										#print('  - product {}: {}'.format(k, Chem.MolToSmiles(product, isomericSmiles = True)))
										try:
											pass
											#product = Chem.MolFromSmiles(Chem.MolToSmiles(product, isomericSmiles = True))
										except Exception as e:
											print('warning: could not draw {}: {}'.format(Chem.MolToSmiles(product, isomericSmiles = USE_STEREOCHEMISTRY), e))

								# Remove untransformed products
								product_set_split = mol_list_to_inchi(outcome).split(' ++ ')
								for reactant in reactants:
									reactant_inchi = mol_list_to_inchi([reactant])
									if reactant_inchi in product_set_split:
										product_set_split.remove(reactant_inchi)
								product_set = ' ++ '.join(product_set_split)

								if product_set not in unique_product_sets:
									unique_product_sets.append(product_set)

						if v: 
							print('\nExpctd {}'.format(product_inchis))
							for product_set in unique_product_sets:
								print('Found  {}'.format(product_set))

						if product_inchis in unique_product_sets:
							if v: print('\nSuccessfully found true products!')
							correct = True
							total_correct += 1
							if len(unique_product_sets) > 1:
								if v: print('...but also found {} more'.format(len(unique_product_sets) - 1))
								# print('Unspecific prediction!')
								# print('ID: {}'.format(example_doc['_id']))
								# print('Reaction: {}'.format(example_doc['RXN_SMILES']))
								# print('Template: {}'.format(rxn_string))
							else:
								total_precise += 1
						else:
							if v: print('\nDid not find true products')
							if len(unique_product_sets) > 1:
								if v: print('...but found {} unexpected ones'.format(len(unique_product_sets)))
						
					if correct or DO_NOT_TEST:	
						# Valid reaction!
						rxn_canonical = canonicalize_transform(rxn_string)
						# print('Pre-resplit: {}'.format(rxn_canonical))
						# Change from inter-molecular to whatever-molecular
						rxn_canonical_split = rxn_canonical.split('>>')
						rxn_canonical = rxn_canonical_split[0][1:-1].replace(').(', '.') + \
							'>>' + rxn_canonical_split[1][1:-1].replace(').(', '.')
						# print('Post resplit: {}'.format(rxn_canonical))
						
						# # TEMPORARY - SKIP INSERTING
						# total_attempted += 1
						# continue

						# Is this old?
						doc = TRANSFORM_DB.find_one({'reaction_smarts': rxn_canonical})

						if doc:
							# Old
							TRANSFORM_DB.update_one(
								{'_id': doc['_id']},
								{
									'$inc': {
										'count': 1,
									},
									'$addToSet': {
										'references': example_doc['_id'],
									},
								}
							)
						else:
							# New
							result = TRANSFORM_DB.insert_one(
								{
									'reaction_smarts': rxn_canonical,
									'rxn_example': reaction_smiles,
									'references': [example_doc['_id']],
									'explicit_H': False,
									'count': 1,
								}
							)
							#print('Created database entry {}'.format(result.inserted_id))
					else:
						print('Did not find true product')
						print('ID: {}'.format(example_doc['_id']))
						print('REACTION: {}'.format(example_doc['RXN_SMILES']))
						print('TEMPLATE: {}'.format(rxn_string))

					total_attempted += 1
			
			except Exception as e:
				print(e)
				if v: 
					print('skipping')
					#raw_input('Enter anything to continue')
				continue


			# Report progress
			if (i % 100) == 0:
				print('{}/{}'.format(i, N))

			# Pause
			#if v: raw_input('Enter anything to continue...')

	except KeyboardInterrupt:
		print('Stopped early!')		
	# except Exception as e:
	# 	print(e)

	total_examples = i + 1
	print('...finished looking through {} reaction records'.format(min([N, i + 1])))

	print('Unmapped reactions in {}/{} ({}%) cases'.format(total_unmapped,
		total_examples, total_unmapped * 100.0 / total_examples))
	print('Partially-mapped reactions in {}/{} ({}%) cases'.format(total_partialmapped,
		total_examples, total_partialmapped * 100.0 / total_examples))
	print('Non-reactions (stereochem only?) in {}/{} ({}%) cases'.format(total_nonreaction,
		total_examples, total_nonreaction * 100.0 / total_examples))
	print('Error-free parsing in {}/{} ({}%) cases'.format(total_attempted, 
		total_examples, total_attempted * 100.0 / total_examples))
	print('Correct product predictions in {}/{} ({}%) cases'.format(total_correct, 
		total_examples, total_correct * 100.0 / total_examples))
	print('Specific product predictions in {}/{} ({}%) cases'.format(total_precise, 
		total_examples, total_precise * 100.0 / total_examples))
	print('Created {} database entries'.format(TRANSFORM_DB.find().count()))

	return True


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing (incl. saving images); defaults to False')
	parser.add_argument('-n', '--num', type = int, default = 50,
						help = 'Maximum number of records to examine; defaults to 50')
	parser.add_argument('-g', '--general', type = str, default = 'n',
						help = 'Use incredibly general templates; defaults to n\n' + 
						'Only appropriate for forward enumeration')
	parser.add_argument('-t', '--test', type = str, default = 'y',
						help = 'Whether to *skip* testing; defaults to y')
	parser.add_argument('--create_mass', type = str, default = 'y',
						help = 'Whether to allow reactions to product mass, defaults y')
	args = parser.parse_args()

	v = args.v
	lg = RDLogger.logger()
	if not v: lg.setLevel(4)
	SUPER_GENERAL_TEMPLATES = args.general in ['y', '1', 'True','Yes', 't', 'T']

	print('SUPER GENERAL TEMPLATES? {}'.format(SUPER_GENERAL_TEMPLATES))
	DO_NOT_TEST = args.test in ['y', '1', 'True','Yes', 't', 'T']
	print('DO NOT TEST? {}'.format(DO_NOT_TEST))
	ALLOW_CREATE_MASS = args.create_mass in ['y', '1', 'True','Yes', 't', 'T']
	print('ALLOW CREATE MASS? {}'.format(ALLOW_CREATE_MASS))

	clear = raw_input('Do you want to clear the {} existing templates? '.format(TRANSFORM_DB.find().count()))
	if clear in ['y', 'Y', 'yes', '1', 'Yes']:
		result = TRANSFORM_DB.delete_many({})
		print('Cleared {} entries from collection'.format(result.deleted_count))

	main(N = args.num)
