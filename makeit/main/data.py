from makeit.utils.neural_fp import *
import numpy as np
import os
import csv
import json


def get_data_full(data_label, config = {}):
	'''This is a helper script to read the data .json object and return
	the training and test data sets separately. This is to allow for an
	already-trained model to be evaluated using the test data (i.e., which
	we know it hasn't seen before)'''

	# Roots
	data_label = data_label.lower()
	data_froot = os.path.join('makeit', 'data')

	###################################################################################
	### WHICH DATASET ARE WE TRYING TO USE?
	###################################################################################

	# Bradley melting point data
	if data_label in ['bradley', 'bradleymp', 'bradley_mp', 'Bradley']:
		dset = 'bradley'
		data_fpath = os.path.join(data_froot, 'BradleyMeltingPointDataset_filtered.csv')
		ftype = 'csv'
		smiles_index = 0
		y_index = 1
		def y_func(x): return x
		y_label = 'Tm (deg C)'

	# Opensol solubility database
	elif data_label in ['opensol', 'open sol']:
		dset = 'opensol'
		data_fpath = os.path.join(data_froot, '20150430SolubilitiesSum_red.csv')
		ftype = 'csv'
		smiles_index = 1
		y_index = 4
		def y_func(x): return np.log10(x)
		# Get solvent
		try:
			temp = config['TRAINING']['solvent']
			solvent = temp
		except:
			solvent = '1-octanol'
		y_label = 'log10({} sol (M))'.format(solvent)

	# Delaney solubility
	elif data_label in ['delaney', 'delaney sol']:
		dset = 'delaney'
		data_fpath = os.path.join(data_froot, 'Delaney2004.txt')
		ftype = 'csv'
		smiles_index = 3
		y_index = 1
		def y_func(x): return x
		y_label = 'log10(aq sol (M))'

	# Ames mutagenicity
	elif data_label in ['ames']:
		dset = 'ames'
		data_fpath = os.path.join(data_froot, 'cas_4337.sdf')
		ftype = 'sdf'
		y_label = 'Ames mutagenicity'

	# Abraham octanol set
	elif data_label in ['abraham', 'abraham sol', 'abraham_oct']:
		dset = 'abraham_oct'
		data_fpath = os.path.join(data_froot, 'AbrahamAcree2014_Octsol_partialSmiles.csv')
		ftype = 'csv'
		smiles_index = 1
		y_index = 5
		def y_func(x): return x
		y_label = 'log10(octanol sol (M))'

	# Other?
	else:
		print('Unrecognized data_label {}'.format(data_label))
		quit(1)


	###################################################################################
	### READ AND TRUNCATE DATA
	###################################################################################

	print('reading data...')
	if ftype == 'csv':
		# Load data from json file
		data = []
		with open(data_fpath, 'r') as data_fid:
			reader = csv.reader(data_fid, delimiter = ',', quotechar = '"')
			for row in reader:
				data.append(row)
	elif ftype == 'sdf':
		data = Chem.SDMolSupplier(data_fpath)
		data_indeces = range(len(data))
	print('done')
		
	# Truncate if necessary
	try:
		truncate_to = int(config['TRAINING']['truncate_to'])
		if ftype == 'csv':
			data = data[:truncate_to]
		elif ftype == 'sdf':
			data_indeces = data_indeces[:truncate_to]
		print('truncated data to first {} samples'.format(truncate_to))
	except:
		pass

	# Get new training_ratio if possible
	try:
		temp = float(config['TRAINING']['ratio'])
		training_ratio = temp
	except:
		training_ratio = 0.8
		pass

	# Get new shuffle seed if possible
	try:
		temp = int(config['TRAINING']['shuffle_seed'])
		np.random.seed(temp)
	except:
		pass
	
	###################################################################################
	### ITERATE THROUGH DATASET AND CREATE NECESSARY DATA LISTS
	###################################################################################

	smiles = []
	mols = []
	y = []
	print('processing data...',)
	# Randomize
	if ftype == 'csv':
		np.random.shuffle(data)
		for i, row in enumerate(data):
			if (i % 1000) == 999:
				print('  {}/{}'.format(i + 1, len(data)))

			try:
				# Filters
				if dset == 'opensol':
					# Only look for ones with the right solvent
					if not (row[2] == solvent):
						continue
					# Do not use entries with DONOTUSE
					if 'DONOTUSE' in row[0]:
						continue
					# Do not use entries with non-positive solubility
					if float(row[y_index]) <= 0:
						continue

				# Molecule first (most likely to fail)
				mol_tensor = molToGraph(Chem.MolFromSmiles(row[smiles_index])).dump_as_tensor()
				this_y = y_func(float(row[y_index]))
				mols.append(mol_tensor)
				y.append(this_y) # Measured log(solubility M/L)
				smiles.append(row[smiles_index]) # Smiles
			except:
				print('Failed to generate graph for {}, y: {}'.format(row[smiles_index], row[y_index]))

	elif ftype == 'sdf':
		np.random.shuffle(data_indeces)
		for j, i in enumerate(data_indeces):
			if (j % 1000) == 999:
				print('  {}/{}'.format(j + 1, len(data_indeces)))
			
			try:
				# Molecule first (most likely to fail)
				mol = suppl[i]
				if not mol:
					print('Failed to generate mol for entry {}'.format(i))
					continue
				mols.append(molToGraph(mol).dump_as_tensor())
				
				if dset == 'ames':
					nonmutagenic = 'nonmutagen' in suppl.GetItemText(i)
					y.append(1 - int(nonmutagenic)) # 1 if mutagenic

				smiles.append(Chem.MolToSmiles(mol)) # Smiles
			except:
				print('Failed to generate mol for entry {}'.format(i))


	###################################################################################
	### PAD MOLECULAR TENSORS
	###################################################################################

	if int(config['TRAINING']['batch_size']) > 1: # NEED TO PAD
		num_atoms = [x.shape[0] for x in mols]
		max_num_atoms = max(num_atoms)
		print('padding tensors up to N_atoms = {}...'.format(max_num_atoms + 1))
		mols = [padGraphTensor(x, max_num_atoms + 1) for x in mols]
	print('done')

	###################################################################################
	### DIVIDE DATA VIA RATIO OR CV
	###################################################################################

	if 'ratio' in config['TRAINING']['data_split']: # split train/notrain
		print('Using first fraction ({}) as training'.format(training_ratio))
		# Create training/development split
		division = int(len(mols) * training_ratio)
		mols_train = mols[:division]
		mols_notrain  = mols[division:]
		y_train = y[:division]
		y_notrain  = y[division:]
		smiles_train = smiles[:division]
		smiles_notrain = smiles[division:]

		# Split notrain up
		mols_val    = mols_notrain[:(len(mols_notrain) / 2)] # first half
		y_val       = y_notrain[:(len(mols_notrain) / 2)] # first half
		smiles_val  = smiles_notrain[:(len(mols_notrain) / 2)] # first half
		mols_test   = mols_notrain[(len(mols_notrain) / 2):] # second half
		y_test      = y_notrain[(len(mols_notrain) / 2):] # second half
		smiles_test = smiles_notrain[(len(mols_notrain) / 2):] # second half
		print('Training size: {}'.format(len(mols_train)))
		print('Validation size: {}'.format(len(mols_val)))
		print('Testing size: {}'.format(len(mols_test)))

	elif 'cv' in config['TRAINING']['data_split']: # cross-validation
		# Default to first fold of 5-fold cross-validation
		folds = 5
		this_fold = 0
		# Read from config file
		try:
			folds = int(config['TRAINING']['folds'].split('/')[1])
			this_fold = int(config['TRAINING']['folds'].split('/')[0]) - 1
		except:
			pass

		# Get target size of each fold
		N = len(mols)
		target_fold_size = int(np.ceil(float(N) / folds))
		# Split up data
		folded_mols 	= [mols[x:x+target_fold_size]   for x in range(0, N, target_fold_size)]
		folded_y 		= [y[x:x+target_fold_size]      for x in range(0, N, target_fold_size)]
		folded_smiles 	= [smiles[x:x+target_fold_size] for x in range(0, N, target_fold_size)]
		print('Split data into {} folds'.format(folds))
		print('...using fold {}'.format(this_fold + 1))

		# Recombine into validation (this_fold), training (everything else), and testing (same as validation)
		mols_train   = [x for fold in (folded_mols[:this_fold] + folded_mols[(this_fold + 1):])     for x in fold]
		y_train      = [x for fold in (folded_y[:this_fold] + folded_y[(this_fold + 1):])           for x in fold]
		smiles_train = [x for fold in (folded_smiles[:this_fold] + folded_smiles[(this_fold + 1):]) for x in fold]
		# Validation is just this fold
		mols_val     = folded_mols[this_fold]
		y_val        = folded_y[this_fold]
		smiles_val   = folded_smiles[this_fold]
		# Test is a copy of validation
		mols_test    = mols_val[:]
		y_test       = y_val[:]
		smiles_test  = smiles_val[:]

	else:
		print('Must specify a data_split type of "ratio" or "cv"')
		quit(1)

	###################################################################################
	### REPACKAGE AS DICTIONARY
	###################################################################################
	
	train = {}; train['mols'] = mols_train; train['y'] = y_train; train['smiles'] = smiles_train
	val   = {}; val['mols']   = mols_val;   val['y']   = y_val;   val['smiles']   = smiles_val
	test  = {}; test['mols']  = mols_test;  test['y']  = y_test;  test['smiles']  = smiles_test

	return (train, val, test)