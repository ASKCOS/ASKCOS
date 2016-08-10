from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import os

def edits_by_count(collection, out = None, minimum = 0.1):
	'''
	Plots the number of edits versus number of examples
	'''
	counts = np.zeros((collection.find().count(), 1)).flatten()
	edits = np.zeros((collection.find().count(), 1)).flatten()
	for (i, doc) in enumerate(collection.find()):
		if 'num_edits' not in doc: continue
		counts[i] = doc['count'] if 'count' in doc else 0
		edits[i] = doc['num_edits'] if 'num_edits' in doc else 0

	print('Total counts: {}'.format(sum(counts)))
	print('Total number of templates: {}'.format(sum(counts != 0)))
	edits = edits[counts >= minimum]
	counts = counts[counts >= minimum]
	print('Total counts using minimum of {}: {}'.format(minimum, sum(counts)))
	counts = np.ones_like(counts) / counts

	
	# Scatterplot
	fig = plt.figure()
	ax = plt.gca()
	ax.scatter(counts, edits, alpha = 0.5)
	# ax.set_yscale('log')
	ax.set_xscale('log')
	ax.axis([np.power(10, np.floor(np.log10(min(counts)))), 1, \
		     0, \
		     100])
	plt.xlabel('1 / template count')
	plt.ylabel('Number of edits')
	plt.title('Transform templates (general_v3) extracted from Lowe database')
	plt.grid(True)
	if out:
		fig.savefig(out + ' edits_count.png')
		np.savetxt(out + ' num_edits gte{}.txt'.format(minimum), zip(counts, edits))
	# plt.show()

	# Histogram
	plt.clf()
	n, bins, patches = plt.hist(edits, range(100), facecolor = 'blue', alpha = 0.5)
	plt.xlabel('Number of edits')
	plt.ylabel('Count')
	plt.title('Histogram of number of edits per template (count >= {})'.format(minimum))
	ymin, ymax = plt.ylim()
	xmin, xmax = plt.xlim()
	plt.axis([0, 100, 0, ymax])
	plt.grid(True)
	plt.savefig(out + ' edit_hist gte{}.png'.format(minimum), bbox_inches = 'tight')

	return

if __name__ == '__main__':
	out_folder = os.path.join(os.path.dirname(__file__), 'output')

	# DATABASE
	from pymongo import MongoClient    # mongodb plugin
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['askcos_transforms']
	# collection_name = 'lowe'
	# collection_name = 'lowe_refs_general'
	# collection_name = 'lowe_refs_general_v2'
	collection_name = 'lowe_refs_general_v3'
	# db = client['reaxys']
	# collection_name = 'transforms_retro_v1'

	collection = db[collection_name]

	edits_by_count(collection, out = os.path.join(out_folder, collection_name), minimum = 50)
