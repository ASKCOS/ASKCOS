from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import os

def get_counts(collection):
	'''
	Gets the 'count' field for all entries in a MongoDB collection
	'''
	counts = np.zeros((collection.find().count(), 1)).flatten()
	for (i, doc) in enumerate(collection.find()):
		counts[i] = doc['count'] if 'count' in doc else 0
	print('Total counts: {}'.format(sum(counts)))
	print('Total number of templates: {}'.format(sum(counts != 0)))
	return counts[counts != 0]

def get_counts_templates(collection):
	'''
	Gets the 'counts' and 'reaction_smarts' field for all entries 
	in a MongoDB collection
	'''
	docs = range(collection.find().count())
	for (i, doc) in enumerate(collection.find()):
		docs[i] = (doc['count'] if 'count' in doc else 0, \
			       doc['reaction_smarts'] if 'reaction_smarts' in doc else '')
	return [x for x in docs if x[0] != 0 and x[1] != '']

def probability_v_rank(counts, out = None):
	'''
	Plots the probability (normalized frequency) versus rank for an
	arbitrary 1D vector of counts
	'''
	counts = np.sort(counts)
	probs = counts / np.sum(counts)
	ranks = range(1, len(probs) + 1)
	ranks.reverse()
	
	# Probability
	fig = plt.figure()
	ax = plt.gca()
	ax.scatter(ranks, probs, alpha = 0.5)
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.axis([1, np.power(10, np.ceil(np.log10(max(ranks)))), \
		     np.power(10, np.floor(np.log10(min(probs)))), \
		     np.power(10, np.ceil(np.log10(max(probs))))])
	plt.xlabel('Rank')
	plt.ylabel('Normalized frequency')
	plt.title('Transform templates extracted from Lowe database')
	plt.grid(True)
	if out:
		fig.savefig(out + ' prob_rank.png')
		np.savetxt(out + ' probs.txt', sorted(probs, reverse = True))
	# plt.show()

	# Count
	fig = plt.figure()
	ax = plt.gca()
	ax.scatter(ranks, counts, alpha = 0.5)
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.axis([1, np.power(10, np.ceil(np.log10(max(ranks)))), \
		     np.power(10, np.floor(np.log10(min(counts)))), \
		     np.power(10, np.ceil(np.log10(max(counts))))])
	plt.xlabel('Rank')
	plt.ylabel('Counts')
	plt.title('Transform templates extracted from Lowe database')
	plt.grid(True)
	if out:
		fig.savefig(out + ' count_rank.png')
		np.savetxt(out + ' counts.txt', sorted(counts, reverse = True))
	# plt.show()

	return

def top_templates(docs, out, n = 10):
	'''
	Finds the top ten transformation templates
	'''
	sorted_docs = sorted(docs, key = lambda x: x[0], reverse = True)
	with open(out + ' top_{}.txt'.format(n), 'w') as fid:
		for i in range(n):
			fid.write('{}\t{}\n'.format(sorted_docs[i][0], sorted_docs[i][1]))
	return 

if __name__ == '__main__':
	out_folder = os.path.join(os.path.dirname(__file__), 'output')

	# DATABASE
	from pymongo import MongoClient    # mongodb plugin
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['askcos_transforms']
	collection_name = 'lowe'
	collection_name = 'lowe_refs_general'
	collection_name = 'lowe_refs_general_v2'
	collection_name = 'lowe_refs_general_v3'
	db = client['reaxys']
	collection_name = 'transforms_retro_v1'

	collection = db[collection_name]

	probability_v_rank(get_counts(collection), out = os.path.join(out_folder, collection_name))
	top_templates(get_counts_templates(collection), os.path.join(out_folder, collection_name), n = 10)
