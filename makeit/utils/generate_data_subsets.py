from makeit.utils.chemnet_connect import * # mongodb connection, gets 'chemicals', 'reactions'
from makeit.utils.database import get_all_ids
from numpy.random import shuffle # for random selection
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
	print '...read chemical IDs'
	shuffle(chemical_ids) 
	print '...shuffled chemical IDs'

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
		if (j % 100) == 0:
			print '{}/{}'.format(j, N)

	print '...constructed data list'

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
	print '...saved json file'

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
	print '...read reaction IDs'
	shuffle(reaction_ids) 
	print '...shuffled reaction IDs'

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
			elif not chemical['mol_weight']: # missing mol_weight
				skip_this_reaction = True
				break
			chemical_names.append(chemical['name'][0])
		if skip_this_reaction:
			continue

		# Append to data
		data.append(chemical_names + [reaction['yield']])
		j = j + 1

		# Report progress
		if (j % 100) == 0:
			print '{}/{}'.format(j, N)

	print '...constructed data list'

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
	print '...saved json file'

	return True

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print('Usage: {} "data type" [max # records]'.format(sys.argv[0]))
		print('    Available [data type]s:')
		print('    - "chemical_names_with_mws"')
		print('    - "reactions_2reac_1prod"')
		quit(1)

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

	else:
		print ('Invalid data type "{}", see usage'.format(sys.argv[1]))
		quit(1)