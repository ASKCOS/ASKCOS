## Meant to use with Daniel Lowe's patent database
# reaction SMILES strings only (for now)

from __future__ import print_function
from numpy.random import shuffle # for random selection
import rdkit.Chem as Chem          # molecule building
from rdkit.Chem import AllChem
import datetime # for info files
import json # for dumping
import sys  # for commanad line
import os   # for file paths

def mols_from_smiles_list(all_smiles):
	'''Given a list of smiles strings, this function creates rdkit
	molecules'''
	mols = []
	for smiles in all_smiles:
		if not smiles: # empty string
			continue
		mols.append(Chem.MolFromSmiles(smiles))
	return mols

def all_bonds_in_mol_list(all_smiles):
	'''Given a list of SMILES strings, this function creates rdkit
	molecules, calculates all bonds present for any molecule in the 
	list, and returns a list of bond tuples, where each tuple 
	characterizes that bond.'''
	bonds = []
	for smiles in all_smiles:
		if not smiles: # empty string
			continue
		mol = Chem.MolFromSmiles(smiles)
		if not mol: # couldn't parse
			print('could not parse {}'.format(smiles))
			continue
		for bond in mol.GetBonds():
			bonds.append(bond_to_label(bond))
	return bonds

def compare_bond_lists(bonds1, bonds2):
	'''This function takes two lists of bond labels and determines the bonds
	that were created and destroyed to get from bonds1 to bonds2'''
	# Look for edge cases
	if (not bonds1) and (not bonds2):
		return ([], [])
	if bonds1 and (not bonds2):
		return ([], bonds1[:])
	if (not bonds1) and bonds2:
		return(bonds2[:], [])

	# Duplicate lists to prevent aliasing issues
	created_bonds = bonds2[:]
	destroyed_bonds = bonds1[:]
	# Remove overlap
	for bond in created_bonds[:]:
		if bond in destroyed_bonds:
			created_bonds.remove(bond)
			destroyed_bonds.remove(bond)
	for bond in destroyed_bonds[:]:
		if bond in created_bonds:
			created_bonds.remove(bond)
			destroyed_bonds.remove(bond)

	return (created_bonds, destroyed_bonds)

def bond_to_label(bond):
	'''This function takes an RDKit bond and creates a label describing
	the most important attributes'''
	atoms = sorted([atom_to_label(bond.GetBeginAtom()), \
				    atom_to_label(bond.GetEndAtom())])

	return atoms[0] + bond.GetSmarts() + atoms[1] 
	return [atoms[0], \
	        atoms[1], \
	        bond.GetSmarts()]

def atom_to_label(atom):
	'''This function takes an RDKit atom and creates a label describing
	the most important attributes'''
	return atom.GetSmarts()
	return '{}[{}]'.format(atom.GetAtomicNum(), \
		                   atom.GetImplicitValence())

def bond_comparison_to_smarts(created_bonds, destroyed_bonds, agents = None):
	'''This function creates a simple reaction SMARTS rule based on the
	created and destroyed bonds in the molcule'''
	# Join components separately
	reactants = '.'.join(destroyed_bonds)
	products  = '.'.join(created_bonds)
	if agents:
		agents = '.'.join(agents)
	else:
		agents = ''
	# Merge 
	return '{}>{}>{}'.format(reactants, agents, products)


def main(db_fpath, N = 15):
	'''Read reactions from Lowe's patent reaction SMILES'''

	# Open file
	data_fid = open(db_fpath, 'r')

	# Look for entries
	for i, line in enumerate(data_fid):

		# Are we done?
		if i == N:
			break

		# Unpack
		reaction_smiles = line.split('\t')[0]
		reaction_smiles = reaction_smiles.split(' ')[0]
		reactants = reaction_smiles.split('>')[0]
		agents    = reaction_smiles.split('>')[1]
		products  = reaction_smiles.split('>')[2]

		# Get bonds
		reactant_bonds = all_bonds_in_mol_list(reactants.split('.'))
		product_bonds  = all_bonds_in_mol_list(products.split('.'))

		# Find how reaction changed bonds
		(created_bonds, destroyed_bonds) = \
			compare_bond_lists(reactant_bonds, product_bonds)

		# Convert to reaction SMARTS (without agent for now)
		rxn_smarts = bond_comparison_to_smarts(created_bonds, destroyed_bonds)

		# Load into rdkit
		rxn = AllChem.ReactionFromSmarts(rxn_smarts)

		# Print
		print('reaction {}'.format(i))
		print('  reaction:  {}'.format(reaction_smiles))
		print('  created:   {}'.format(created_bonds))
		if (len(destroyed_bonds) - len(created_bonds)) > 3:
			print('  destroyed: ** missing product, many bonds lost')
			continue
		else:
			print('  destroyed: {}'.format(destroyed_bonds))

		# Test reaction
		print('reaction test')
		print('  {}'.format(rxn_smarts))
		reactants = mols_from_smiles_list(reactants.split('.'))
		products = mols_from_smiles_list(products.split('.'))
		outcomes = rxn.RunReactants(reactants)
		for j, outcome in outcomes:
			print('- outcome {}/{}'.format(j, len(outcomes)))
			for k, product in outcome:
				print('  product {}: {}'.format(k, Chem.MolToSmiles(product)))
		
		# Report progress
		if (i % 1000) == 0:
			print('{}/{}'.format(i, N))

	print('...finished looking through data')

	return True


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print('Usage: {} "data.rsmi" [max # records]'.format(sys.argv[0]))
		quit(1)

	# Run 
	if len(sys.argv) == 3:
		main(sys.argv[1], N = int(sys.argv[2]))
	else:
		main(sys.argv[1])