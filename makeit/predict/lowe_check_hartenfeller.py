# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import argparse
import numpy as np     	      	   # for simple calculations
import os                          # for saving
import sys
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from pymongo import MongoClient    # mongodb plugin
from bson import ObjectId
from makeit.predict.summarize_reaction_outcome import summarize_reaction_outcome
from makeit.embedding.descriptors import edits_to_vectors, oneHotVector # for testing
import re
import time
from tqdm import tqdm
import cPickle as pickle
import itertools 


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing; defaults to False')
	parser.add_argument('--candidate_collection', type = str, default = 'candidate_edits_8_9_16', 
						help = 'Collection of candidates to write to; defaults to candidate_edits_8_9_16')
	parser.add_argument('--singleonly', type = bool, default = True,
						help = 'Whether to record major product only; defaults to True')
	parser.add_argument('--check', type = bool, default = False,
						help = 'Whether to check already-done ones or not')
	
	args = parser.parse_args()
	v = bool(args.v)

	from rdkit import RDLogger
	lg = RDLogger.logger()
	lg.setLevel(4)

	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['prediction']
	candidates = db[args.candidate_collection]

	done_ids = set()
	if args.check:
		for doc in candidates.find({'found_hartenfeller': {'$exists': True}}, []):
			done_ids.add(doc['_id'])
		print('Checked completed entries')
		print('{} records have hartenfeller checked'.format(len(done_ids)))

	
	fid_templates = open('makeit/predict/hartenfeller.txt', 'r')
	templates = []
	for line in fid_templates.readlines():
		line = line.replace('"', '')
		rxn = AllChem.ReactionFromSmarts(line.split('\t')[1])
		templates.append(rxn)

	print('Loaded {} templates'.format(len(templates)))

	fid_data = open('makeit/predict/lowe_data_edits/lowe_data.pickle', 'rb')
	fid_labels = open('makeit/predict/lowe_data_edits/lowe_labels.pickle', 'rb')


	legend_labels = pickle.load(fid_labels) # first doc is legend
	CANDIDATE_SMILES = legend_labels['candidate_smiles']
	CANDIDATE_EDITS  = legend_labels['candidate_edits_compact']
	REACTION_TRUE    = legend_labels['reaction_true']
	RXDID            = legend_labels['rxdid']

	legend_data = pickle.load(fid_data) # first doc is legend
	CANDIDATE_EDITS_COMPACT = legend_data['candidate_edits_compact']
	ATOM_DESC_DICT          = legend_data['atom_desc_dict']
	REACTION_TRUE_ONEHOT    = legend_data['reaction_true_onehot']

	N_samples =  legend_data['N_examples']

	for i in tqdm(range(N_samples)):
		doc_data = pickle.load(fid_data)
		doc_label = pickle.load(fid_labels)

		if ObjectId(doc_label[RXDID]) in done_ids: continue

		rxn_smiles = '['.join(doc_label[REACTION_TRUE].split('[')[:-1])

		reactant_smiles = rxn_smiles.split('>')[0]
		product_smiles  = rxn_smiles.split('>')[2]
		
		# Load reactant molecules
		reactants = [Chem.MolFromSmiles(x) for x in reactant_smiles.split('.')]

		found_true = False
		for template in templates:
			if found_true: break

			# Try all combinations of reactants that fit template
			combinations = itertools.combinations(reactants, rxn.GetNumReactantTemplates())
			unique_product_sets = []
			for combination in combinations:

				# Perform transformation
				try:
					outcomes = template.RunReactants(list(combination))
				except Exception as e:
					continue
				if not outcomes: continue # no match
				for j, outcome in enumerate(outcomes):
					outcome = outcome[0] # all products represented as single mol by transforms

					# Reduce to largest (longest) product only?
					candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY)

					candidate_smiles = max(candidate_smiles.split('.'), key = len)
					outcome = Chem.MolFromSmiles(candidate_smiles)
					if not outcome: continue

					# Remove mapping before matching
					[x.ClearProp('molAtomMapNumber') for x in outcome.GetAtoms() if x.HasProp('molAtomMapNumber')] # remove atom mapping from outcome
					if Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY) == product_smiles:
						if v: print('Matched true [{}]'.format(product_smiles))
						found_true = True


		# Prepare doc and insert
		if v:
			if found_true: 
				print('Found true product')
			else:
				print('Did not find true product')

		try:

			res = candidates.update_one(
				{'_id': ObjectId(doc_label[RXDID])},
				{'$set':
					{'found_hartenfeller': found_true}
				}
			)
		except Exception as e:
			print(e)
			continue

	total_processed = candidates.count({'found_hartenfeller': {'$exists': True}})
	total_found     = candidates.count({'found_hartenfeller': True})

	print('Out of {} processed, {} were found by Hartenfeller templates'.format(total_processed, total_found))