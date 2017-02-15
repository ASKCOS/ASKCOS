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
from collections import defaultdict

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
	padUpTo = len(candidate_list)
	N_e1 = 20
	N_e2 = 20
	N_e3 = 20
	N_e4 = 20

	N_e1_trim = 1
	N_e2_trim = 1
	N_e3_trim = 1
	N_e4_trim = 1

	# Initialize
	x_h_lost = np.zeros((1, padUpTo, N_e1, F_atom))
	x_h_gain = np.zeros((1, padUpTo, N_e2, F_atom))
	x_bond_lost = np.zeros((1, padUpTo, N_e3, F_bond))
	x_bond_gain = np.zeros((1, padUpTo, N_e4, F_bond))

	# Get reactant descriptors
	atom_desc_dict = edits_to_vectors([], reactants, return_atom_desc_dict = True)

	# Populate arrays
	for (c, edits) in enumerate(candidate_edits):
		if any([len(edit) > N_e for edit in edits]):
			continue
		if c == padUpTo: break
		edit_h_lost_vec, edit_h_gain_vec, \
			edit_bond_lost_vec, edit_bond_gain_vec = edits_to_vectors(edits, reactants, atom_desc_dict = atom_desc_dict)

		N_e1_trim = max(N_e1_trim, len(edit_h_lost_vec))
		N_e2_trim = max(N_e2_trim, len(edit_h_gain_vec))
		N_e3_trim = max(N_e3_trim, len(edit_bond_lost_vec))
		N_e4_trim = max(N_e4_trim, len(edit_bond_gain_vec))

		for (e, edit_h_lost) in enumerate(edit_h_lost_vec):
			if e >= N_e1: continue
			x_h_lost[0, c, e, :] = edit_h_lost
		for (e, edit_h_gain) in enumerate(edit_h_gain_vec):
			if e >= N_e2: continue
			x_h_gain[0, c, e, :] = edit_h_gain
		for (e, edit_bond_lost) in enumerate(edit_bond_lost_vec):
			if e >= N_e3: continue
			x_bond_lost[0, c, e, :] = edit_bond_lost
		for (e, edit_bond_gain) in enumerate(edit_bond_gain_vec):
			if e >= N_e4: continue
			x_bond_gain[0, c, e, :] = edit_bond_gain

	# Trim down
	x_h_lost = x_h_lost[:, :, :N_e1_trim, :]
	x_h_gain = x_h_gain[:, :, :N_e2_trim, :]
	x_bond_lost = x_bond_lost[:, :, :N_e3_trim, :]
	x_bond_gain = x_bond_gain[:, :, :N_e4_trim, :]

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

	pred = model.predict(xs + [xc[:, 7:], xc[:, 1:7], xc[:, 0]], batch_size = 5, verbose = True)[0]

	def softmax(x):
		"""Compute softmax values for each sets of scores in x."""
		e_x = np.exp(x - np.max(x))
		return e_x / e_x.sum()

	pred_scaled = softmax(pred)

	# Turn edit predictions into molecule predictions
	candidate_dict = defaultdict(float)
	for i, p in enumerate(pred_scaled):
		if i >= len(candidate_list): break
		candidate_dict[candidate_list[i][0]] += p



	rank = ss.rankdata(pred)

	fname = raw_input('Enter file name to save to: ')

	# Edit level
	with open(os.path.join(FROOT, fname + ' edits.dat'), 'w') as fid:
		fid.write('FOR REACTANTS {} WITH CONTEXT {}\n'.format(Chem.MolToSmiles(reactants), context))
		fid.write('Candidate product\tCandidate edit\tScore\tProb\tRank\n')
		for (c, candidate) in enumerate(candidate_list):
			candidate_smile = candidate[0]
			candidate_edit = candidate[1]
			fid.write('{}\t{}\t{}\t{}\t{}\n'.format(
				candidate_smile, candidate_edit, pred[c], pred_scaled[c], 1 + len(pred) - rank[c]
			))
	print('Wrote to file {}'.format(os.path.join(FROOT, fname + ' edits.dat')))

	# Molecule-level
	with open(os.path.join(FROOT, fname + ' mols.dat'), 'w') as fid:
		fid.write('FOR REACTANTS {} WITH CONTEXT {}\n'.format(Chem.MolToSmiles(reactants), context))
		fid.write('Candidate product\tProbability\tRank\n')
		for (c, candidate) in enumerate(sorted(candidate_dict.iteritems(), reverse = True, key = lambda x: x[1])):
			candidate_smile = candidate[0]
			candidate_prob = candidate[1]
			fid.write('{}\t{}\t{}\n'.format(
				candidate_smile, candidate_prob, c + 1
			))
	print('Wrote to file {}'.format(os.path.join(FROOT, fname + ' mols.dat')))

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('--tag', type = str,
						help = 'Tag for model to load from')
	parser.add_argument('--mincount', type = int, default = 100,
						help = 'Mincount of templates, default 100')
	parser.add_argument('--Nc', type = int, default = 1000,
						help = 'Number of candidates to consider, default 1000')
	parser.add_argument('--Nh1', type = int, default = 200,
						help = 'Number of hidden nodes in first layer, default 200')
	parser.add_argument('--Nh2', type = int, default = 100,
						help = 'Number of hidden nodes in second layer, default 100')
	parser.add_argument('--Nh3', type = int, default = 50,
						help = 'Number of hidden nodes in third layer, ' + 
								'immediately before summing, default 50')
	parser.add_argument('--Nhf', type = int, default = 50,
						help = 'Number of hidden nodes in layer between summing ' +
								'and final score, default 50')
	parser.add_argument('--l2', type = float, default = 0.0,
						help = 'l2 regularization parameter for each Dense layer, default 0.0')
	parser.add_argument('--lr', type = float, default = 0.01, 
						help = 'Learning rate, default 0.01')
	parser.add_argument('--context_weight', type = float, default = 100.0,
						help = 'Weight assigned to contextual effects, default 100.0')
	parser.add_argument('--enhancement_weight', type = float, default = 0.1,
					help = 'Weight assigned to enhancement factor, default 0.1')
	args = parser.parse_args()


	# Labels
	tag = args.tag
	FROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'output')
	FROOT = os.path.join(FROOT, tag)
	if not os.path.isdir(FROOT):
		os.mkdir(FROOT)
	MODEL_FPATH = os.path.join(FROOT, 'model.json')
	WEIGHTS_FPATH = os.path.join(FROOT, 'weights.h5')
	HIST_FPATH = os.path.join(FROOT, 'hist.csv')
	TEST_FPATH = os.path.join(FROOT, 'probs.dat')
	HISTOGRAM_FPATH = os.path.join(FROOT, 'histogram %s.png')
	ARGS_FPATH = os.path.join(FROOT, 'args.json')

	with open(ARGS_FPATH, 'r') as fid:
		print('This args: {}'.format(args.__dict__))
		import json
		old_dict = json.load(fid)
		for (key, val) in old_dict.iteritems():
			args.__dict__[key] = val
		print('Reloaded args: {}'.format(args.__dict__))

	mol = Chem.MolFromSmiles('[C:1][C:2]')
	(a, _, b, _) = edits_to_vectors((['1'],[],[('1','2',1.0)],[]), mol)
	F_atom = len(a[0])
	F_bond = len(b[0])

	template_collection = 'transforms_forward_v1'
	mincount = int(args.mincount)
	singleonly = True
	v = False
	N_h1 = int(args.Nh1)
	N_h2 = int(args.Nh2)
	N_h3 = int(args.Nh3)
	N_hf = int(args.Nhf)
	l2v = float(args.l2)
	lr = float(args.lr)
	N_c = int(args.Nc) # number of candidate edit sets - UNUSED
	N_e = 5 # maximum number of edits per class - UNUSED
	context_weight = float(args.context_weight)
	enhancement_weight = float(args.enhancement_weight)
	optimizer          = args.optimizer
	inner_act          = args.inner_act
	TARGET_YIELD       = bool(args.yd)

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
	Transformer.load(templates, mincount = mincount, get_retro = False, get_synth = True, lowe = False)
	print('Out of {} database templates,'.format(templates.count()))
	print('Loaded {} templates'.format(Transformer.num_templates))
	Transformer.reorder()
	print('Sorted by count, descending')

	# Load models
	from makeit.predict.reaxys_score_candidates_from_edits_compact_v2 import build
	model = build(F_atom = F_atom, F_bond = F_bond, N_h1 = N_h1, 
				N_h2 = N_h2, N_h3 = N_h3, N_hf = N_hf, l2v = l2v, inner_act = inner_act,
				context_weight = context_weight, enhancement_weight = enhancement_weight, TARGET_YIELD = TARGET_YIELD,
				absolute_score = True)

	model.load_weights(WEIGHTS_FPATH, by_name = True)
	print('Loaded weights from file')

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
				print('Try one of the following: {}'.format(', '.join([doc['name'] for doc in SOLVENT_DB.find() if 'name' in doc])))
				continue
			solvent_vec = [doc['c'], doc['e'], doc['s'], doc['a'], doc['b'], doc['v']]
			solvent = doc['name']

			# Unreacting reagents
			reagents = [Chem.MolFromSmiles(reagent) for reagent in raw_input('Enter SMILES of reagents: ').split('.')]
			if None in reagents:
				print('Could not parse all reagents!')
				continue
			if not reagents:
				continue
			reagent_fp = np.zeros(256)
			for reagent in reagents:
				reagent_fp += np.array(AllChem.GetMorganFingerprintAsBitVect(reagent, 2, nBits = 256))

			# Combine
			xc = np.resize(np.array([T] + solvent_vec + list(reagent_fp)), (1, 263))

			context = 'solv:{}, rgt:{}, T:{}'.format(
				solvent, '.'.join([Chem.MolToSmiles(reagent) for reagent in reagents]), T
			)

			print('Number of reactant atoms: {}'.format(len(reactants.GetAtoms())))
			# Report current reactant SMILES string
			[a.ClearProp('molAtomMapNumber') for a in reactants.GetAtoms() if a.HasProp('molAtomMapNumber')]
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