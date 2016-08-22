from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
# matplotlib.rcParams.update({'font.size': 22})
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
from collections import defaultdict
import copy
import os

if __name__ == '__main__':
	out_folder = os.path.join(os.path.dirname(__file__), 'output')

	# DATABASE
	from pymongo import MongoClient    # mongodb plugin
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaxys']
	REDUCTANT_DB = db['reductants']
	
	reductants = defaultdict(int) # dict of counts
	for doc in REDUCTANT_DB.find({}, ['_id', 'SMILES', 'count']):
		reductants[doc['SMILES']] = doc['count']

	counts = reductants.values()
	labels = reductants.keys()
	total_count = sum(counts)

	zipsorted = sorted(zip(counts, labels), key = lambda x: x[0])
	counts = np.array([x[0] for x in zipsorted])
	labels = [x[1] for x in zipsorted]
	probs = counts * 1.0 / np.sum(counts)
	ranks = range(1, len(probs) + 1)
	ranks.reverse()

	# Oxidizing agent scatterplot
	fig = plt.figure(figsize=(6,4), dpi = 300)
	ax = plt.gca()
	ax.scatter(ranks, probs, alpha = 0.8)
	ax.set_yscale('log')
	ax.set_xscale('log')
	plt.xlabel('Rank')
	plt.ylabel('Probability')
	ax.axis([1, np.power(10, np.ceil(np.log10(max(ranks)))), \
	     np.power(10, np.floor(np.log10(min(probs)))), \
	     np.power(10, np.ceil(np.log10(max(probs))))])
	plt.title('Suspected reductants in Reaxys')
	plt.grid(True)
	N_top = 5
	for (i, (label, x, y)) in enumerate(zip(labels, ranks, probs)[-N_top:]):
		plt.annotate(
			label, 
			xy = (x, y), xytext = (1, 1 + (i - N_top + 0.5) / N_top),
			textcoords = 'axes fraction', ha = 'right', va = 'center',
			bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.8),
			arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
	fig.savefig(os.path.join(out_folder, 'reductants.png'))
