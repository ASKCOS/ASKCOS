from __future__ import print_function
from makeit.utils.chemnet_connect import * # mongodb connection, gets 'chemicals', 'reactions'
from makeit.utils.database import get_all_ids
from numpy.random import shuffle # for random selection
import rdkit.Chem as Chem          # molecule building
import rdkit.Chem.FunctionalGroups as rdFGs
import datetime # for info files
import json # for dumping
import sys  # for commanad line
import os   # for file paths

def get_data_folder_root():
	'''Returns folder path where all data is stored'''
	return os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')

def dump_to_data_file(data, label, details = ''):
	'''This function acts as a save function for the other generation functions
	in the module; it uses json.dump to save data sets. An info file is 
	also created which contains, at minimum, a timestamp.'''

	# Build path names
	data_folder = get_data_folder_root()
	data_fname = os.path.join(data_folder, '{}.json'.format(label))
	info_fname = os.path.join(data_folder, '{}.info'.format(label))

	# Write to data file
	data_fid = open(data_fname, 'w')
	json.dump(data, data_fid)
	data_fid.close()

	# Write to info file
	info_fid = open(info_fname, 'w')
	time_now = datetime.datetime.utcnow()
	info_fid.write('{}.json generated at UTC {}\n\n'.format(label, time_now))
	info_fid.write('File details\n------------\n')
	info_fid.write(details)
	info_fid.close()

	return True

def chemical_names_with_mws(N = 50000):
	'''Sample chemicals collection for chemicals with a valid name and valid
	molecular weight. Currently, this function does not attempt to fill in 
	missing database information.

	Saved data is [str(name), float(mol_weight)]'''
	
	# Randomize list of chemical IDs
	chemical_ids = get_all_ids(chemicals)
	print('...read chemical IDs')
	shuffle(chemical_ids) 
	print('...shuffled chemical IDs')

	# Look for valid entries
	data = []
	j = 0 # successful entries == len(data)
	for chemical_id in chemical_ids:

		# Are we done?
		if j == N:
			break

		# Find entry
		chemical = chemicals.find_one({'_id' : chemical_id})

		# Filter
		if not chemical['name']: # missing name
			continue
		if not chemical['mol_weight']: # missing mol_weight
			continue

		# Append to list
		data.append([chemical['name'][0], chemical['mol_weight']])
		j = j + 1

		# Report progress
		if (j % 1000) == 0:
			print('{}/{}'.format(j, N))

	print('...constructed data list')

	# Write details
	details = 'Found {} random chemicals from database'.format(len(data))
	details += ' satisfying the following criteria:\n'
	details += '- chemical[\'name\'] == True\n'
	details += '- chemical[\'mol_weight\'] == True\n'
	details += '\nData list consists of entries:\n'
	details += '  [str(name), float(mol_weight)]\n'

	# Save
	dump_to_data_file(data, 'chemical_names_with_mws_{}'.format(len(data)), 
		details = details)
	print('...saved json file')

	return True

def reactions_2reac_1prod(N = 10000):
	'''Sample database for reactions with exactly two reactants and one 
	product (A + B -> C). Chemicals are checked to make sure they have 
	a non-empty name field. 

	Saved data is [str(A_name), str(B_name), str(C_name), float(yield)]'''

	# Build filter (2 reactants and 1 product)
	db_filter = {'reactants' : {'$size' : 2}, 'products' : {'$size' : 1}}

	# Randomize list of chemical IDs
	reaction_ids = get_all_ids(reactions, db_filter = db_filter)
	print('...read reaction IDs')
	shuffle(reaction_ids) 
	print('...shuffled reaction IDs')

	# Look for valid entries
	data = []
	j = 0 # successful entries == len(data)
	for reaction_id in reaction_ids:

		# Are we done?
		if j == N:
			break

		# Find entry
		reaction = reactions.find_one({'_id' : reaction_id})

		# Unpack
		[idA, idB] = reaction['reactants']
		[idC]      = reaction['products']
		chemical_names = []

		# Filter
		if not reaction['yield']: # missing yield
			continue
		skip_this_reaction = False
		for chemical_id in [idA, idB, idC]:
			chemical = chemicals.find_one({'_id' : chemical_id})
			if not chemical: # ID not found
				skip_this_reaction = True
				break
			elif not chemical['name']: # missing name
				skip_this_reaction = True
				break
			chemical_names.append(chemical['name'][0])
		if skip_this_reaction:
			continue

		# Append to data
		data.append(chemical_names + [reaction['yield']])
		j = j + 1

		# Report progress
		if (j % 1000) == 0:
			print('{}/{}'.format(j, N))

	print('...constructed data list')

	# Write details
	details = 'Found {} random reactions from database'.format(len(data))
	details += ' satisfying the following criteria:\n'
	details += '- exactly 2 reactants\n'
	details += '- exactly 1 product\n'
	details += 'for the associated chemicals...\n'
	details += '- chemical[\'name\'] == True\n'
	details += '\nData list consists of entries:\n'
	details += '  [str(A name), str(B name), str(C name), float(yield)]\n'
	
	# Save
	dump_to_data_file(data, 'reactions_2reac_1prod_{}'.format(len(data)), 
		details = details)
	print('...saved json file')

	return True


def reactions_2reac_rdsmiles(N = 10000):
	'''Sample database for reactions with exactly two reactants. Chemicals
	are checked to make sure they can be parsed by RDKit by SMILES. Yields
	of exactly 0.75 are filtered out

	Saved data is [str(SMILES), str(SMILES), float(yield)]'''

	# Build filter (2 reactants and 1 product)
	db_filter = {'reactants' : {'$size' : 2}, 'yield' : {'$ne' : 0.75}}

	# Randomize list of reaction IDsgit gui
	reaction_ids = get_all_ids(reactions, db_filter = db_filter)
	print('...read reaction IDs')
	shuffle(reaction_ids) 
	print('...shuffled reaction IDs')

	# Look for valid entries
	data = []
	j = 0 # successful entries == len(data)
	for reaction_id in reaction_ids:

		# Are we done?
		if j == N:
			break

		# Find entry
		reaction = reactions.find_one({'_id' : reaction_id})

		# Unpack
		[idA, idB] = reaction['reactants']
		chemical_names = []

		# Filter
		if not reaction['yield']: # missing yield
			continue
		skip_this_reaction = False
		for chemical_id in [idA, idB]:
			chemical = chemicals.find_one({'_id' : chemical_id})
			if not chemical: # ID not found
				skip_this_reaction = True
				break
			elif not chemical['smiles']: # missing smiles
				skip_this_reaction = True
				break
			this_mol = Chem.MolFromSmiles(chemical['smiles'])
			if not this_mol:
				skip_this_reaction = True
				break
			chemical_names.append(Chem.MolToSmiles(this_mol))
		if skip_this_reaction:
			continue

		# Append to data
		data.append(chemical_names + [reaction['yield']])
		j = j + 1

		# Report progress
		if (j % 1000) == 0:
			print('{}/{}'.format(j, N))

	print('...constructed data list')

	# Write details
	details = 'Found {} random reactions from database'.format(len(data))
	details += ' satisfying the following criteria:\n'
	details += '- exactly 2 reactants\n'
	details += '- yield != 0.75 exactly\n'
	details += 'for the associated chemicals...\n'
	details += '- chemical[\'smiles\'] == True\n'
	details += '- can be parsed by RDKit'
	details += '\nData list consists of entries:\n'
	details += '  [str(SMILES), str(SMILES), float(yield)]\n'
	
	# Save
	dump_to_data_file(data, 'reactions_2reac_rdsmiles_{}'.format(len(data)), 
		details = details)
	print('...saved json file')

	return True

def reactions_2reac_smilesyield(N = 10000):
	'''Sample database for reactions with exactly two reactants. Yields of 
	exactly 0.75 are filtered out.

	Saved data is [str(RXN_SMILES), float(yield)]'''

	# Build filter (2 reactants and 1 product)
	db_filter = {'yield' : {'$ne' : 0.75}}

	# Randomize list of reaction IDsgit gui
	reaction_ids = get_all_ids(reactions, db_filter = db_filter)
	print('...read reaction IDs')
	shuffle(reaction_ids) 
	print('...shuffled reaction IDs')

	# Look for valid entries
	data = []
	j = 0 # successful entries == len(data)
	for reaction_id in reaction_ids:

		# Are we done?
		if j == N:
			break

		# Find entry
		reaction = reactions.find_one({'_id' : reaction_id})

		# Filter
		if not reaction['yield']: # missing yield
			continue
		if not reaction['smiles']: # missing reaction
			continue

		# Append to data
		data.append([reaction['smiles'], reaction['yield']])
		j = j + 1

		# Report progress
		if (j % 1000) == 0:
			print('{}/{}'.format(j, N))

	print('...constructed data list')

	# Write details
	details = 'Found {} random reactions from database'.format(len(data))
	details += ' satisfying the following criteria:\n'
	details += '- yield != 0.75 exactly\n'
	details += '\nData list consists of entries:\n'
	details += '  [str(RXN_SMILES), float(yield)]\n'
	
	# Save
	dump_to_data_file(data, 'reactions_2reac_smilesyield_{}'.format(len(data)), 
		details = details)
	print('...saved json file')

	return True

def chemical_names(N = 100000):
	'''Sample chemicals collection for chemicals with a valid name, meant to
	generate a dataset suitable for fitting a tokenizer.

	Saved data is str(name)'''
	
	# Randomize list of chemical IDs
	chemical_ids = get_all_ids(chemicals)
	print('...read chemical IDs')
	shuffle(chemical_ids) 
	print('...shuffled chemical IDs')

	# Look for valid entries
	data = []
	j = 0 # successful entries == len(data)
	for chemical_id in chemical_ids:

		# Are we done?
		if j == N:
			break

		# Find entry
		chemical = chemicals.find_one({'_id' : chemical_id})

		# Filter
		if not chemical['name']: # missing name
			continue

		# Append to list
		data.append(chemical['name'][0])
		j = j + 1

		# Report progress
		if (j % 5000) == 0:
			print('{}/{}'.format(j, N))

	print('...constructed data list')

	# Write details
	details = 'Found {} random chemicals from database'.format(len(data))
	details += ' satisfying the following criteria:\n'
	details += '- chemical[\'name\'] == True\n'
	details += '\nData list consists of entries:\n'
	details += '  str(name)\n'

	# Save
	dump_to_data_file(data, 'chemical_names_{}'.format(len(data)), 
		details = details)
	print('...saved json file')

	return True
	
def chemical_rdsmiles(N = 100000):
	'''Sample chemicals collection for chemicals that can be parsed into RDKit

	Saved data is str(SMILES)'''
	
	# Randomize list of chemical IDs
	chemical_ids = get_all_ids(chemicals)
	print('...read chemical IDs')
	shuffle(chemical_ids) 
	print('...shuffled chemical IDs')

	# Look for valid entries
	data = []
	j = 0 # successful entries == len(data)
	for chemical_id in chemical_ids:

		# Are we done?
		if j == N:
			break

		# Find entry
		chemical = chemicals.find_one({'_id' : chemical_id})

		# Try loading smiles
		if chemical['smiles']:
			try:
				this_mol = Chem.MolFromSmiles(chemical['smiles'])
				if not this_mol:
					continue
			except:
				continue

		# Append to list
		data.append(Chem.MolToSmiles(this_mol))
		j = j + 1

		# Report progress
		if (j % 5000) == 0:
			print('{}/{}'.format(j, N))

	print('...constructed data list')

	# Write details
	details = 'Found {} random chemicals from database'.format(len(data))
	details += ' satisfying the following criteria:\n'
	details += '- could be parsed by RDKit using SMILES\n'
	details += '\nData list consists of entries:\n'
	details += '  str(SMILES)\n'

	# Save
	print(data)
	dump_to_data_file(data, 'chemical_rdmols_{}'.format(len(data)), 
		details = details)
	print('...saved json file')

	return True

def chemical_rdsmiles_rdfuncgroups(N = 10000):
	'''Sample chemicals collection for chemicals that can be parsed into RDKit

	Saved data is str(SMILES), list(funcgroup FP)'''
	

	# Get hierarchy
	hierarchy = rdFGs.BuildFuncGroupHierarchy()

	# Randomize list of chemical IDs
	chemical_ids = get_all_ids(chemicals)
	print('...read chemical IDs')
	shuffle(chemical_ids) 
	print('...shuffled chemical IDs')

	# Look for valid entries
	data = []
	j = 0 # successful entries == len(data)
	for chemical_id in chemical_ids:

		# Are we done?
		if j == N:
			break

		# Find entry
		chemical = chemicals.find_one({'_id' : chemical_id})

		# Try loading smiles
		if chemical['smiles']:
			try:
				this_mol = Chem.MolFromSmiles(chemical['smiles'])
				this_fp = rdFGs.CreateMolFingerprint(this_mol, hierarchy)
				if not this_mol:
					continue
			except:
				continue

		# Append to list
		data.append([Chem.MolToSmiles(this_mol), this_fp])
		j = j + 1

		# Report progress
		if (j % 100) == 0:
			print('{}/{}'.format(j, N))

	print('...constructed data list')

	# Write details
	details = 'Found {} random chemicals from database'.format(len(data))
	details += ' satisfying the following criteria:\n'
	details += '- could be parsed by RDKit using SMILES\n'
	details += '\nData list consists of entries:\n'
	details += '  str(SMILES), list(FunctionalGroup FP)\n'

	# Save
	dump_to_data_file(data, 'chemical_rdsmiles_rdfuncgroups_{}'.format(len(data)), 
		details = details)
	print('...saved json file')

	return True

def all_reaction_dois():
	'''This function searches for all DOIs references in the reaction
	datababse.'''

	# Look for valid entries
	data = []
	N = reactions.find().count()
	num_refs = 0
	num_rxns_with_refs = 0
	num_rxns_with_doi = 0
	for i, reaction in enumerate(reactions.find()):
		this_data = []

		# Get rxid tag
		rxid = reaction['rxid'][0]
		
		# Try reading reference list
		try: 
			references = reaction['rx' + str(rxid)]
			num_rxns_with_refs += 1
			# Get list of dois from all references
			if references:
				for reference in references:
					num_refs += 1
					try:
						this_data += reference['doi']
					except:
						pass
		except:
			pass

		# Append to data
		if this_data:
			num_rxns_with_doi += 1
			data += this_data

		# Report progress
		if (i % 10000) == 0:
			print('{}/{}'.format(i, N))

	print('...constructed data list ({} DOIs from {} refs)'.format(len(data), num_refs))

	# Write details
	details = 'Found {} DOIs in database '.format(len(data))
	details += ' from {} references'.format(num_refs)
	details += ', representing {} different reactions'.format(num_rxns_with_refs)
	details += ', but only {} of those had DOI information'.format(num_rxns_with_doi)

	# Save
	dump_to_data_file(data, 'all_reaction_dois_{}'.format(len(data)), 
		details = details)
	print('...saved json file')

	return True

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print('Usage: {} "data type" [max # records]'.format(sys.argv[0]))
		print('    Available [data type]s:')
		print('    - "chemical_names_with_mws"')
		print('    - "reactions_2reac_1prod"')
		print('    - "chemical_names"')
		print('    - "all_reaction_dois"')
		print('    - "chemical_rdsmiles"')
		print('    - "reactions_2reac_smilesyield"')
		print('    - "chemical_rdsmiles_rdfuncgroups"')

	# Molecular weight training set
	if sys.argv[1] == 'chemical_names_with_mws':
		if len(sys.argv) == 3:
			chemical_names_with_mws(int(sys.argv[2]))
		else:
			chemical_names_with_mws()

	elif sys.argv[1] == 'reactions_2reac_1prod':
		if len(sys.argv) == 3:
			reactions_2reac_1prod(int(sys.argv[2]))
		else:
			reactions_2reac_1prod()

	elif sys.argv[1] == 'reactions_2reac_rdsmiles':
		if len(sys.argv) == 3:
			reactions_2reac_rdsmiles(int(sys.argv[2]))
		else:
			reactions_2reac_rdsmiles()

	elif sys.argv[1] == 'chemical_names':
		if len(sys.argv) == 3:
			chemical_names(int(sys.argv[2]))
		else:
			chemical_names()

	elif sys.argv[1] == 'all_reaction_dois':
		all_reaction_dois()

	elif sys.argv[1] == 'chemical_rdsmiles':
		if len(sys.argv) == 3:
			chemical_rdsmiles(int(sys.argv[2]))
		else:
			chemical_rdsmiles()		

	elif sys.argv[1] == 'reactions_2reac_smilesyield':
		if len(sys.argv) == 3:
			reactions_2reac_smilesyield(int(sys.argv[2]))
		else:
			reactions_2reac_smilesyield()		

	elif sys.argv[1] == 'chemical_rdsmiles_rdfuncgroups':
		if len(sys.argv) == 3:
			chemical_rdsmiles_rdfuncgroups(int(sys.argv[2]))
		else:
			chemical_rdsmiles_rdfuncgroups()

	else:
		print('Invalid data type "{}", see usage'.format(sys.argv[1]))
		quit(1)