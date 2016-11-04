# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import numpy as np
import os
import sys
import argparse
import h5py # needed for save_weights, fails otherwise
from keras import backend as K 
from keras.models import Sequential, Model, model_from_json
from keras.optimizers import *
import makeit.retro.transformer as transformer 
from makeit.retro.canonicalization import SmilesFixer
from pymongo import MongoClient    # mongodb plugin
from makeit.predict.summarize_reaction_outcome import summarize_reaction_outcome
from makeit.embedding.descriptors import edits_to_vectors, oneHotVector # for testing
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import theano.tensor as T
import cPickle as pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt    # for visualization
import scipy.stats as ss
import itertools
import time
from tqdm import tqdm

def reactants_to_candidate_edits(reactants):
	candidate_list = [] # list of tuples of (product smiles, edits required)
	for template in tqdm(Transformer.templates):
		# Perform transformation
		try:
			outcomes = template['rxn_f'].RunReactants([reactants])
		except Exception as e:
			if v: print(e)
			continue
		if not outcomes: continue # no match
		for j, outcome in enumerate(outcomes):
			outcome = outcome[0] # all products represented as single mol by transforms
			try:
				outcome.UpdatePropertyCache()
				Chem.SanitizeMol(outcome)
				[a.SetProp('molAtomMapNumber', a.GetProp('old_molAtomMapNumber')) \
					for (i, a) in enumerate(outcome.GetAtoms()) \
					if 'old_molAtomMapNumber' in a.GetPropsAsDict()]
			except Exception as e:
				if v: print(e)
				continue
			if v: print('Outcome SMILES: {}'.format(Chem.MolToSmiles(outcome)))

			
			
			# Reduce to largest (longest) product only?
			candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY)
			if singleonly: 
				candidate_smiles = max(candidate_smiles.split('.'), key = len)
				outcome = Chem.MolFromSmiles(candidate_smiles)
				
			# Find what edits were made
			edits = summarize_reaction_outcome(reactants, outcome)
			if v: print(edits)

			# Remove mapping before matching
			[x.ClearProp('molAtomMapNumber') for x in outcome.GetAtoms() if x.HasProp('molAtomMapNumber')] # remove atom mapping from outcome

			# Overwrite candidate_smiles without atom mapping numbers
			candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY)

			# Add to ongoing list
			if (candidate_smiles, edits) not in candidate_list:
				candidate_list.append((candidate_smiles, edits))

	return candidate_list 

def preprocess_candidate_edits(reactants, candidate_list):
	candidate_smiles = [a for (a, b) in candidate_list]
	candidate_edits =  [b for (a, b) in candidate_list]

	print('Generated {} unique edit sets'.format(len(candidate_list)))

	# Initialize
	x_h_lost = np.zeros((1, padUpTo, N_e, F_atom))
	x_h_gain = np.zeros((1, padUpTo, N_e, F_atom))
	x_bond_lost = np.zeros((1, padUpTo, N_e, F_bond))
	x_bond_gain = np.zeros((1, padUpTo, N_e, F_bond))

	# Populate arrays
	for (c, edits) in enumerate(candidate_edits):
		if any([len(edit) > 5 for edit in edits]):
			continue
		if c == padUpTo: break
		edit_h_lost_vec, edit_h_gain_vec, \
			edit_bond_lost_vec, edit_bond_gain_vec = edits_to_vectors(edits, reactants)
		for (e, edit_h_lost) in enumerate(edit_h_lost_vec):
			x_h_lost[0, c, e, :] = edit_h_lost
		for (e, edit_h_gain) in enumerate(edit_h_gain_vec):
			x_h_gain[0, c, e, :] = edit_h_gain
		for (e, edit_bond_lost) in enumerate(edit_bond_lost_vec):
			x_bond_lost[0, c, e, :] = edit_bond_lost
		for (e, edit_bond_gain) in enumerate(edit_bond_gain_vec):
			x_bond_gain[0, c, e, :] = edit_bond_gain

	# Get rid of NaNs
	x_h_lost[np.isnan(x_h_lost)] = 0.0
	x_h_gain[np.isnan(x_h_gain)] = 0.0
	x_bond_lost[np.isnan(x_bond_lost)] = 0.0
	x_bond_gain[np.isnan(x_bond_gain)] = 0.0
	x_h_lost[np.isinf(x_h_lost)] = 0.0
	x_h_gain[np.isinf(x_h_gain)] = 0.0
	x_bond_lost[np.isinf(x_bond_lost)] = 0.0
	x_bond_gain[np.isinf(x_bond_gain)] = 0.0

	return [x_h_lost, x_h_gain, x_bond_lost, x_bond_gain]

def score_candidates(reactants, candidate_list, xs, context, xc):

	pred = model.predict(xs + [xc], batch_size = 20)[0]
	rank = ss.rankdata(pred)

	fname = raw_input('Enter file name to save to: ') + '.dat'
	with open(os.path.join(FROOT, fname), 'w') as fid:
		fid.write('FOR REACTANTS {} WITH CONTEXT {}\n'.format(Chem.MolToSmiles(reactants, context)))
		fid.write('Candidate product\tCandidate edit\tProbability\tRank\n')
		for (c, candidate) in enumerate(candidate_list):
			if c >= padUpTo: break
			candidate_smile = candidate[0]
			candidate_edit = candidate[1]
			fid.write('{}\t{}\t{}\t{}\n'.format(
				candidate_smile, candidate_edit, pred[c], 1 + len(pred) - rank[c]
			))
	print('Wrote to file {}'.format(os.path.join(FROOT, fname)))

if __name__ == '__main__':
	FROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'output')

	parser = argparse.ArgumentParser()
	parser.add_argument('--tag', type = str,
						help = 'Tag for model to load from')
	parser.add_argument('--mincount', type = int, default = 50,
						help = 'Mincount of templates, default 50')

	args = parser.parse_args()


	mol = Chem.MolFromSmiles('[C:1][C:2]')
	(a, _, b, _) = edits_to_vectors((['1'],[],[('1','2',1.0)],[]), mol)
	F_atom = len(a[0])
	F_bond = len(b[0])

	tag = args.tag
	template_collection = 'transforms_forward_v1'
	mincount = int(args.mincount)
	singleonly = True
	padUpTo = 100
	N_e = 5
	v = False

	# Silence warnings
	from rdkit import RDLogger
	lg = RDLogger.logger()
	lg.setLevel(4)

	# Load transformer
	from pymongo import MongoClient    # mongodb plugin
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaxys']
	templates = db[template_collection]
	SOLVENT_DB = db['solvents']
	Transformer = transformer.Transformer()
	Transformer.load(templates, mincount = mincount, get_retro = False, get_synth = True, lowe = True)
	print('Out of {} database templates,'.format(templates.count()))
	print('Loaded {} templates'.format(Transformer.num_templates))
	Transformer.reorder()
	print('Sorted by count, descending')

	# Load models
	model = model_from_json(open(os.path.join(FROOT, 'model{}.json'.format(tag))).read())
	model.compile(loss = 'categorical_crossentropy',
			optimizer = SGD(lr = 0.0),
			metrics = ['accuracy']
	)
	model.load_weights(os.path.join(FROOT, 'weights{}.h5'.format(tag)))
	print('Loaded model from file')

	# Wait for user prompt
	while True:
		try:
			print('-------------------')
			# 1) Reactants
			reactants = Chem.MolFromSmiles(raw_input('Enter SMILES of reactants: '))
			if not reactants:
				print('Could not parse!')
				continue
			# 2) Temperature
			T = float(raw_input('Enter temperature [C]: '))

			# 3) Solvent
			solvent = raw_input('Enter solvent (SMILES or name): ')
			solvent_mol = Chem.MolFromSmiles(solvent)
			if solvent_mol: 
				doc = SOLVENT_DB.find_one({'_id': Chem.MolToSmiles(solvent_mol)})
			else:
				doc = SOLVENT_DB.find_one({'name': solvent})
			if not doc:
				print('Could not parse solvent!')
				print('Try one of the following: {}'.format(', '.join([doc['name'] for doc in SOLVENT_DB.find()])))
				continue
			solvent_vec = [doc['c'], doc['e'], doc['s'], doc['a'], doc['b'], doc['v']]
			solvent = doc['name']

			# Unreacting reagents
			reagents = [Chem.MolFromSmiles(reagent) for reagent in raw_input('Enter SMILES of reagents: ').split('.')]
			if None in reagents:
				print('Could not parse all reagents!')
				continue
			reagent_fp = np.zeros(256)
			for reagent in reagents:
				reagent_fp += np.array(AllChem.GetMorganFingerprintAsBitVect(reagent, 2, nBits = 256))

			# Combine
			xc = np.array([T] + solvent_vec + list(reagent_fp))
			context = 'solv:{}, rgt:{}, T:{}'.format(
				solvent, '.'.join([Chem.MolToSmiles(reagent) for reagent in reagents]), T
			)

			print('Number of reactant atoms: {}'.format(len(reactants.GetAtoms())))
			# Report current reactant SMILES string
			print('Reactants w/o map: {}'.format(Chem.MolToSmiles(reactants)))
			# Add new atom map numbers
			[a.SetProp('molAtomMapNumber', str(i+1)) for (i, a) in enumerate(reactants.GetAtoms())]
			# Report new reactant SMILES string
			print('Reactants w/ map: {}'.format(Chem.MolToSmiles(reactants)))

			# Generate candidates
			candidate_list = reactants_to_candidate_edits(reactants)
			# Convert to matrices
			xs = preprocess_candidate_edits(reactants, candidate_list)
			# Score and save to file
			score_candidates(reactants, candidate_list, xs, context, xc)
		except Exception as e:
			print(e)
			print('')