## Meant to use with Daniel Lowe's patent database
# reaction SMILES strings only (for now)

from __future__ import print_function
from numpy.random import shuffle # for random selection
import rdkit.Chem as Chem          # molecule building
from rdkit.Chem import AllChem
from collections import defaultdict
import rdkit.Chem.Draw as Draw
import datetime # for info files
import json # for dumping
import sys  # for commanad line
import os   # for file paths
import re 
import itertools

def mols_from_smiles_list(all_smiles):
	'''Given a list of smiles strings, this function creates rdkit
	molecules'''
	mols = []
	for smiles in all_smiles:
		if not smiles: continue
		mols.append(Chem.MolFromSmiles(smiles))
	return mols

# def all_bonds_in_mol_list(all_smiles):
# 	'''Given a list of SMILES strings, this function creates rdkit
# 	molecules, calculates all bonds present for any molecule in the 
# 	list, and returns a list of bond tuples, where each tuple 
# 	characterizes that bond.'''
# 	bonds = []
# 	for smiles in all_smiles:
# 		if not smiles: # empty string
# 			continue
# 		mol = Chem.MolFromSmiles(smiles)
# 		if not mol: # couldn't parse
# 			if v: print('could not parse {}'.format(smiles))
# 			continue
# 		for bond in mol.GetBonds():
# 			bonds.append(bond_to_label(bond))
# 	return bonds

# def compare_bond_lists(bonds1, bonds2):
# 	'''This function takes two lists of bond labels and determines the bonds
# 	that were created and destroyed to get from bonds1 to bonds2'''
# 	# Look for edge cases
# 	if (not bonds1) and (not bonds2):
# 		return ([], [])
# 	if bonds1 and (not bonds2):
# 		return ([], bonds1[:])
# 	if (not bonds1) and bonds2:
# 		return(bonds2[:], [])

# 	# Duplicate lists to prevent aliasing issues
# 	created_bonds = bonds2[:]
# 	destroyed_bonds = bonds1[:]
# 	# Remove overlap
# 	for bond in created_bonds[:]:
# 		if bond in destroyed_bonds:
# 			created_bonds.remove(bond)
# 			destroyed_bonds.remove(bond)
# 	for bond in destroyed_bonds[:]:
# 		if bond in created_bonds:
# 			created_bonds.remove(bond)
# 			destroyed_bonds.remove(bond)

# 	return (created_bonds, destroyed_bonds)

# def bond_to_label(bond):
# 	'''This function takes an RDKit bond and creates a label describing
# 	the most important attributes'''
# 	atoms = sorted([atom_to_label(bond.GetBeginAtom()), \
# 				    atom_to_label(bond.GetEndAtom())])

# 	return atoms[0] + bond.GetSmarts() + atoms[1] 
# 	return [atoms[0], \
# 	        atoms[1], \
# 	        bond.GetSmarts()]

# def atom_to_label(atom):
# 	'''This function takes an RDKit atom and creates a label describing
# 	the most important attributes'''
# 	return atom.GetSmarts()
# 	return '{}[{}]'.format(atom.GetAtomicNum(), \
# 		                   atom.GetImplicitValence())

# def bond_comparison_to_smarts(created_bonds, destroyed_bonds, agents = None):
# 	'''This function creates a simple reaction SMARTS rule based on the
# 	created and destroyed bonds in the molcule'''
# 	# Join components separately


# 	# Look until reaction can be initialized
# 	while True:

# 		# Identify which atom indeces show up multiple times
# 		destroyed_bonds_atoms = [[] for x in destroyed_bonds]
# 		mapped_atoms = []
# 		for i, destroyed_bond in enumerate(destroyed_bonds):
# 			atom_indeces = re.findall(r'\:([0-9]*)\]', destroyed_bond)
# 			destroyed_bonds_atoms[i] += [int(x) for x in set(atom_indeces)]
# 			mapped_atoms += [int(x) for x in set(atom_indeces)]

# 		# Get duplicates
# 		duplicates = set([x for x in mapped_atoms if mapped_atoms.count(x) > 1])
# 		if v: print('duplicates: {}'.format(duplicates))

# 		if not duplicates:
# 			break

# 		# Get the first conflict
# 		atom_index = duplicates.pop()
# 		if v: print('first duplicate: {}'.format(atom_index))
# 		group = '('
# 		bonds_in_group = []
# 		for i, atom_indeces in enumerate(destroyed_bonds_atoms):
# 			if atom_index in atom_indeces:
# 				group += destroyed_bonds[i] + '.'
# 				bonds_in_group.append(i)
# 		group = group[:-1] + ')'

# 		# Update destroyed_bonds to reflect grouping
# 		destroyed_bonds[bonds_in_group[0]] = group
# 		for bond_index in bonds_in_group[1:]:
# 			destroyed_bonds.remove(destroyed_bonds[bond_index])
		
# 	# Combine now
# 	reactants = '.'.join(destroyed_bonds)
# 	products  = '.'.join(created_bonds)
# 	if agents:
# 		agents = '.'.join(agents)
# 	else:
# 		agents = ''

# 	# Merge 
# 	return '{}>{}>{}'.format(reactants, agents, products)

def map_reaction(reactants, products):
	'''Re-assigns atom numbers for reactants and products to reflect most likely
	atom mapping'''
	# TODO
	return reactants, products

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

	if atom1.GetAtomicNum() != atom2.GetAtomicNum(): return True # must be true for atom mapping
	if atom1.GetTotalNumHs() != atom2.GetTotalNumHs(): return True
	if atom1.GetFormalCharge() != atom2.GetFormalCharge(): return True
	if atom1.GetDegree() != atom2.GetDegree(): return True
	# TODO: add # pi electrons like ICSynth?

	if level >= 1:
		# Check neighbors too
		neighbors1 = sorted([atom.GetAtomicNum() for atom in atom1.GetNeighbors()])
		neighbors2 = sorted([atom.GetAtomicNum() for atom in atom2.GetNeighbors()])
		if neighbors1 != neighbors2: return True

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
		err = 1
	if len(prod_atoms) != len(reac_atoms):
		if v: print('warning: total number of tagged atoms differ, stoichometry != 1?')
		err = 1

	# Find differences 
	changed_atoms = []
	changed_atom_tags = []
	# Product atoms that are different from reactant atom equivalent
	for i, prod_tag in enumerate(prod_atom_tags):
		for j, reac_tag in enumerate(reac_atom_tags):
			if reac_tag != prod_tag: continue
			if reac_tag not in changed_atom_tags: # don't bother comparing if we know this atom changes
				if atoms_are_different(prod_atoms[i], reac_atoms[j]):
					changed_atoms.append(reac_atoms[j])
					changed_atom_tags.append(reac_tag)
					break
	# Reactant atoms that do not appear in product (tagged leaving groups)
	for j, reac_tag in enumerate(reac_atom_tags):
		if reac_tag not in changed_atom_tags:
			if reac_tag not in prod_atom_tags:
				changed_atoms.append(reac_atoms[j])
				changed_atom_tags.append(reac_tag)

	if v: print('{} tagged atoms in reactants change 1-atom properties'.format(len(changed_atom_tags)))
	for smarts in [atom.GetSmarts() for atom in changed_atoms]:
		if v: print('  {}'.format(smarts))

	return changed_atoms, changed_atom_tags, err

def draw_reaction_smiles(rxn_smiles, fpath = 'test.png'):

	opts = Draw.DrawingOptions()
	opts.elemDict = defaultdict(lambda: (0,0,0))
	opts.noCarbonSymbols = False
	opts.selectColor = (1, 0, 0)

	rxn = AllChem.ReactionFromSmarts(rxn_smiles)
	img = Draw.ReactionToImage(rxn, subImgSize = (300, 300), highlightAtoms = [], options = opts)
	img.save(fpath)

	return

def expand_atoms_to_use(mol, atoms_to_use):
	'''Given an RDKit molecule and a list of AtomIdX which should be included
	in a fragment, this function expands the list of AtomIdXs to include one 
	nearest neighbor'''

	new_atoms_to_use = atoms_to_use[:]
	for atom in mol.GetAtoms():
		if atom.GetIdx() in atoms_to_use: continue
		for neighbor in atom.GetNeighbors():
			if neighbor.GetIdx() in atoms_to_use:
				new_atoms_to_use.append(atom.GetIdx())
				break
	return new_atoms_to_use

def get_fragments_for_changed_atoms(mols, changed_atom_tags, radius = 0):
	'''Given a list of RDKit mols and a list of changed atom tags, this function
	computes the SMILES string of molecular fragments using MolFragmentToSmiles 
	for all changed fragments'''
	fragments = ''
	for mol in mols:
		atoms_to_use = []
		for atom in mol.GetAtoms():
			# Check self (only tagged atoms)
			if ':' in atom.GetSmarts():
				if atom.GetSmarts().split(':')[1][:-1] in changed_atom_tags:
					atoms_to_use.append(atom.GetIdx())
					continue
		# Check neighbors (any atom)
		for k in range(radius):
			atoms_to_use = expand_atoms_to_use(mol, atoms_to_use)

		if not atoms_to_use: continue
		fragments += '(' + AllChem.MolFragmentToSmiles(mol, atoms_to_use) + ').'
	return fragments[:-1]

def mol_list_to_inchi(mols):
	'''List of RDKit molecules to InChI string separated by ++'''
	inchis = [Chem.MolToInchi(mol) for mol in mols]
	return ' ++ '.join(sorted(inchis))

def mol_list_from_inchi(inchis):
	'''InChI string separated by ++ to list of RDKit molecules'''
	return [Chem.MolFromInchi(inchi.strip()) for inchi in inchis.split('++')]

def main(db_fpath, N = 15, radius = 1):
	'''Read reactions from Lowe's patent reaction SMILES'''

	# Open file
	data_fid = open(db_fpath, 'r')

	# Define scoring variables
	total_templates = 0 # total reactions simulated (excludes skipped)
	total_correct = 0 # actual products predicted
	total_precise = 0 # ONLY actual products predicted

	out_fid = open('test/transforms/all_transforms.txt', 'w')

	# Look for entries
	for i, line in enumerate(data_fid):

		# Are we done?
		if i == N:
			break

		if v: print('##################################')
		if v: print('###        RXN {}'.format(i))
		if v: print('##################################')

		# Unpack
		reaction_smiles = line.split('\t')[0].split(' ')[0] # remove fragment portion, too
		reactants, agents, products = [mols_from_smiles_list(x) for x in 
										[mols.split('.') for mols in reaction_smiles.split('>')]]
		if None in reactants + agents + products: 
			if v: print('Could not parse all molecules in reaction, skipping')
			continue

		# Map atoms
		reactants, products = map_reaction(reactants, products)

		# Calculate changed atoms
		changed_atoms, changed_atom_tags, err = get_changed_atoms(reactants, products)
		if err: 
			if v: print('skipping')
			continue

		# Draw
		if v: print(reaction_smiles)
		#draw_reaction_smiles(reaction_smiles, fpath = 'test/transforms/rxn_{}.png'.format(i))

		# Get fragments
		if v: print('Using radius {} to build fragments'.format(radius))
		reactant_fragments = get_fragments_for_changed_atoms(reactants, changed_atom_tags, radius = radius)
		product_fragments  = get_fragments_for_changed_atoms(products,  changed_atom_tags, radius = radius)

		# Report transform
		rxn_string = '{}>>{}'.format(reactant_fragments, product_fragments)
		if v: print('Overall fragment transform: {}'.format(rxn_string))

		# Load into RDKit
		rxn = AllChem.ReactionFromSmarts(rxn_string)

		# Analyze
		if v: print('Original number of reactants: {}'.format(len(reactants)))
		if v: print('Number of reactants in transform: {}'.format(rxn.GetNumReactantTemplates()))

		# Run reaction
		try:
			# Try all combinations of reactants that fit template
			combinations = itertools.combinations(reactants, rxn.GetNumReactantTemplates())
			unique_product_sets = []
			for combination in combinations:
				outcomes = rxn.RunReactants(list(combination))
				if not outcomes: continue
				if v: print('\nFor reactants {}'.format([Chem.MolToSmiles(mol) for mol in combination]))
				for j, outcome in enumerate(outcomes):
					if v: print('- outcome {}/{}'.format(j + 1, len(outcomes)))
					for k, product in enumerate(outcome):
						if v: print('  product {}: {}'.format(k, Chem.MolToSmiles(product)))
					product_set = mol_list_to_inchi(outcome)
					if product_set not in unique_product_sets:
						unique_product_sets.append(product_set)

			if v: print('\nExpctd {}'.format(mol_list_to_inchi(products)))
			for product_set in unique_product_sets:
				if v: print('Found  {}'.format(product_set))

			if mol_list_to_inchi(products) in unique_product_sets:
				if v: print('\nSuccessfully found true products!')
				total_correct += 1
				if len(unique_product_sets) > 1:
					if v: print('...but also found {} more'.format(len(unique_product_sets) - 1))
				else:
					total_precise += 1
					out_fid.write(rxn_string + '\n')
			else:
				if v: print('\nDid not find true products')
				if len(unique_product_sets) > 1:
					if v: print('...but found {} unexpected ones'.format(len(unique_product_sets)))
			
			total_templates += 1
		
		except Exception as e:
			if v: print(e)
			if v: print('skipping')
			continue


		# Report progress
		if (i % 1000) == 0:
			print('{}/{}'.format(i, N))

		# Pause
		# raw_input('Enter anything to continue...')

	print('...finished looking through {} reaction records'.format(N))

	print('Correct product predictions in {}/{} ({}%) cases'.format(total_correct, 
		total_templates, total_correct * 100.0 / total_templates))
	print('Specific product predictions in {}/{} ({}%) cases'.format(total_precise, 
		total_templates, total_precise * 100.0 / total_templates))

	out_fid.close()

	return True


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print('Usage: {} "data.rsmi" [max # records]'.format(sys.argv[0]))
		quit(1)

	# Verbose?
	if '-v' in sys.argv: 
		v = True
	else:
		v = False

	# Run 
	if len(sys.argv) >= 4:
		main(sys.argv[1], N = int(sys.argv[2]), radius = int(sys.argv[3]))
	elif len(sys.argv) >= 3:
		main(sys.argv[1], N = int(sys.argv[2]))
	else:
		main(sys.argv[1])

