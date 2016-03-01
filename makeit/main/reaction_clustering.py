from __future__ import print_function
from makeit.utils.neural_fp import *
from sklearn.cluster import DBSCAN
import rdkit.Chem as Chem
import numpy as np
import json
import sys
import os


def generate_rxn_embeddings_and_save(embedding_func, save_fpath):
	'''This is a helper script to read the data .json object and return
	the embedded reactions'''


	# Load data from json file
	print('reading data...',)
	data = json.load('makeit/data/reactions_rdsmiles_968173.json')
	print('done')
		
	# Truncate if necessary
	data = data[:5000]

	# Get an example embedding so we know the size
	mol_tensor = molToGraph(Chem.MolFromSmiles(data[0][0][0])).dump_as_tensor()
	print(mol_tensor.shape)
	single_mol_as_array = np.array([mol_tensor])
	embedding = embedding_func([single_mol_as_array])

	# Parse data into individual components
	smiles = []
	embeddings = np.zeros(len(data), len(embedding))
	print('processing data...',)
	solvents = {}
	# Randomize
	for i, row in enumerate(data):
		if (i % 100) == 99:
			print('  {}/{}'.format(i + 1, len(data) - 1))
		try:
			reacstring = ''
			for reactant in row[0]:
				mol_tensor = molToGraph(Chem.MolFromSmiles(reactant)).dump_as_tensor()
				embeddings[i, :] -= embedding_func([np.array([mol_tensor])])
				reacstring += reactant + '.'
			reacstring = reacstring[:-1] + '>>'
			for product in row[1]:
				mol_tensor = molToGraph(Chem.MolFromSmiles(product)).dump_as_tensor()
				embeddings[i, :] += embedding_func([np.array([mol_tensor])])
				reacstring += product + '.'
			reacstring = reacstring[:-1]
			smiles.append(reacstring)
		except:
			print('Failed to generate graph for [reacs, prods] = {}'.format(row))

	json.dump(embeddings, save_fpath + '_embeddings.json')
	json.dump(smiles, save_fpath + '_smiles.json')
	print('saved embeddings and smiles for reactions')

	return (embebddings, smiles)

def cluster(embeddings, smiles, eps = 0.3, min_samples = 10):
	# Create DBSCAN object
	db = DBSCAN(eps = eps, min_samples = min_samples).fit(embeddings)
	labels = db.labels_

	# Number of clusters in labels, ignoring noise if present.
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
	print('Estimated number of clusters: %d' % n_clusters_)

	for label in set(labels):
		cluster_reactions = (labels == label)
		print('Looking at cluster label {} ({} members)'.format(label, sum(cluster_reactions)))
		for i, e in enumerate(smiles):
			if cluster_reactions(i):
				print('{}. {}'.format(i, e))