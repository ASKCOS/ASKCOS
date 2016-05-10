## Meant to use with Daniel Lowe's patent database
# reaction SMILES strings only (for now)

from __future__ import print_function
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
from xml.dom import minidom

# DATABASE
from pymongo import MongoClient    # mongodb plugin
client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaction_examples']
collection = db['lowe_1976-2013_USPTOgrants_reactions']

def doc_to_dic(document):
	'''Converts a parsed XML file to a dictionary using relevant fields'''
	dic = {}
	dic['reaction_smiles'] = \
		str(document.getElementsByTagName('dl:reactionSmiles')[0].firstChild.nodeValue)
	products = []
	for node in document.getElementsByTagName('product'):
		products.append(mol_to_dic(node, withAmounts = True))
	reactants = []
	for node in document.getElementsByTagName('reactant'):
		reactants.append(mol_to_dic(node))
	catalysts = []; solvents = []; spectators = [];
	for node in document.getElementsByTagName('spectator'):
		role = node.attributes.getNamedItem('role').value
		if role == 'catalyst':
			catalysts.append(mol_to_dic(node))
		elif role == 'solvent':
			solvents.append(mol_to_dic(node))
		else:
			raise ValueError('Unknown spectator role: {}'.format(role))
	if not reactants or not products:
		raise ValueError('No reactants/products?')
	dic['reactants'] = reactants
	dic['products'] = products 
	if catalysts: dic['catalysts'] = catalysts 
	if solvents: dic['solvents'] = solvents 
	if spectators: dic['spectators'] = spectators

	# Make yield more accessible
	if sum(['yield' in product for product in products]) == 1:
		dic['yield'] = [product['yield'] for product in products if 'yield' in product][0]

	return dic

def mol_to_dic(node, withAmounts = False):
	'''Converts a node containing molecule information into a
	dictionary'''
	dic = {}
	# Get name
	dic['name'] = str(node.getElementsByTagName('name')[0].firstChild.nodeValue)
	# If exact entity match, more data is available
	#print(node.toprettyxml())
	#entityType = node.getElementsByTagName('dl:entityType')[0].firstChild.nodeValue
	#if entityType == 'exact' or entityType == 'definiteReference':
	identifiers = {
		child.attributes.getNamedItem('dictRef').value : \
		child.attributes.getNamedItem('value').value \
		for child in node.getElementsByTagName('identifier')
	}
	if 'cml:inchi' in identifiers.keys():
		mol = AllChem.MolFromInchi(str(identifiers['cml:inchi']))
	elif 'cml:smiles' in identifiers.keys():
		mol = AllChem.MolFromSmiles(str(identifiers['cml:smiles']))
	else:
		print('identifiers: {}'.format(identifiers.keys()))
		raise ValueError('No molecular identifier for {}'.format(dic['name']))
	if not mol: raise ValueError('Couldnt parse molecule: {}'.format(identifiers))

	Chem.SanitizeMol(mol)
	dic['smiles'] = AllChem.MolToSmiles(mol, isomericSmiles=True)
	dic['inchi'] = AllChem.MolToInchi(mol)
	# elif entityType == 'chemicalClass':
	# 	pass # name is all we get
	# else:
	# 	raise ValueError('Unknown entityType for molecule: {}'.format(entityType))
	# Quantity?
	if withAmounts:
		amounts = {
			child.attributes.getNamedItem('units').value : \
			child.firstChild.nodeValue \
			for child in node.getElementsByTagName('amount')
		}
		if 'unit:percentYield' in amounts.keys():
			dic['yield'] = float(amounts['unit:percentYield'])
		if 'unit:g' in amounts.keys():
			dic['amount(g)'] = float(amounts['unit:g'])
	return dic

def main(db_fpath, N = 15):
	'''Read reactions from Lowe's patent reaction SMILES'''

	try:
		# Open file
		file_generator = get_reaction_file(db_fpath)
		print(file_generator)
		documents = []
		for i, rxn in enumerate(file_generator):
			if i == N:
				break

			print('~~~~~~~ {} ~~~~~~'.format(i))
			print('{}: {}'.format(i, rxn))
			document = minidom.parse(rxn)
			try:
				dic = doc_to_dic(document)
				documents.append(dic)
			except ValueError as e:
				print(e)

			# Report progress and insert every 1000
			if ((i+1) % 1000) == 0:
				print('{}/{}'.format(i+1, N))
				result = collection.insert(documents)
				documents = []

		if documents: result = collection.insert(documents)
	except KeyboardInterrupt:
		print('Stopped early!')		

	print('Created {} database entries'.format(collection.find().count()))

	return True


def get_reaction_file(fpath):
	'''This function is used to recursively traverse a directory and find
	all .cml files contained within. It is used as a generator.'''
	for dir, subdirs, files in os.walk(fpath):
		subdirs.sort(reverse = True)
		for subdir in subdirs:
			for p in get_reaction_file(subdir):
				yield p
		files.sort(reverse = True)
		for file in files:
			yield os.path.join(dir, file)

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('data_file', type = str, 
		 				help = 'File where each line is an atom-mapped smiles reaction')
	parser.add_argument('-n', '--num', type = int, default = 50,
						help = 'Maximum number of records to load; defaults to 50')
	args = parser.parse_args()

	clear = raw_input('Do you want to clear the {} existing examples? '.format(collection.find().count()))
	if clear in ['y', 'Y', 'yes', '1', 'Yes']:
		result = collection.delete_many({})
		print('Cleared {} entries from collection'.format(result.deleted_count))

	main(args.data_file, N = args.num)

