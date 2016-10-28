from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import sys
import os
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import re
import itertools

def atom_to_ox_state(atom):
	known_rules = ['K', 'Na', 'Cl', 'Br', 'F', 'I', 'O', 'H']
	if atom.GetSymbol in ['K', 'Na']:
		return +1
	if atom.GetSymbol() in ['Cl', 'Br', 'F', 'I']:
		# diatom exception
		if any([a.GetSymbol == atom.GetSymbol() for a in atom.GetNeighbors()]):
			return 0
		return -1
	if atom.GetSymbol() == 'O':
		# peroxide exception
		if any([a.GetSymbol() == 'O' for a in atom.GetNeighbors()]): 
			return -1 
		# OF exception
		if any([a.GetSymbol() == 'F' for a in atom.GetNeighbors()]):
			return +2
		return -2 + atom.GetTotalNumHs() + atom.GetFormalCharge()
	elif atom.GetSymbol() == 'H':
		return +1
	else: # mainly C or N
		ox = 0
		for a in atom.GetNeighbors():
			if a.GetSymbol() in known_rules:
				ox -= atom_to_ox_state(a) # offset neighbors
		return ox - atom.GetTotalNumHs() + atom.GetFormalCharge()

def get_oxidation_change(reactant, product):
	if not reactant: raise ValueError('Supplied reactant is NoneType')
	if not product: raise ValueError('Supplied product is NoneType')

	# Get atom maps that appear in product
	maps = []
	for atom in product.GetAtoms():
		if atom.HasProp('molAtomMapNumber'):
			maps.append(atom.GetProp('molAtomMapNumber'))
	
	reactant_ox = 0
	for atom in reactant.GetAtoms():
		if atom.HasProp('molAtomMapNumber'):
			if atom.GetProp('molAtomMapNumber') in maps:
				a_ox = atom_to_ox_state(atom)
				# print('Atom {}, oxidation {}'.format(atom.GetSmarts(), a_ox))
				reactant_ox += a_ox
	# print('TOTAL OF MAPPED ATOMS IN REACTANTS: {}'.format(reactant_ox))

	product_ox = 0
	for atom in product.GetAtoms():
		if atom.HasProp('molAtomMapNumber'):
			if atom.GetProp('molAtomMapNumber') in maps:
				a_ox = atom_to_ox_state(atom)
				# print('Atom {}, oxidation {}'.format(atom.GetSmarts(), a_ox))
				product_ox += a_ox
	# print('TOTAL OF MAPPEPD ATOMS IN PRODUCTS: {}'.format(product_ox))

	delta = float(product_ox - reactant_ox)
	# print('CHANGE: {}'.format(delta))
	# if delta == 0:
	# 	print('NEUTRAL! {}'.format(reaction['RXN_SMILES']))
	# elif delta < 0:
	# 	print('### REDUCTION ### {}'.format(reaction['RXN_SMILES']))
	# elif delta > 0:
	# 	print('### OXIDATION ###  {}'.format(reaction['RXN_SMILES']))
	return delta