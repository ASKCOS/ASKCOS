from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
# matplotlib.rcParams.update({'font.size': 22})
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
from collections import defaultdict
import copy
import os
from tqdm import tqdm

if __name__ == '__main__':
	out_folder = os.path.join(os.path.dirname(__file__), 'output')

	# DATABASE
	from pymongo import MongoClient    # mongodb plugin
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaxys']
	INSTANCE_DB = db['instances']
	CHEMICAL_DB = db['chemicals']
	
	counts = defaultdict(int)
	names = defaultdict(str)
	smiles = defaultdict(str)
	mapped_ids = set()

	for doc in tqdm(INSTANCE_DB.find({}, ['RXD_SOLXRN'], no_cursor_timeout = True)):
		for xrn in doc['RXD_SOLXRN']:
			counts[xrn] += 1
			if xrn not in mapped_ids:
				chem_doc = CHEMICAL_DB.find_one({'_id': xrn})
				if not chem_doc:
					smiles[xrn] = '?? {}'.format(xrn)
					names[xrn] = '?? {}'.format(xrn)
				else:
					smiles[xrn] = chem_doc['SMILES']
					names[xrn] = chem_doc['IDE_CN']
				mapped_ids.add(xrn)
	N_rxdids = INSTANCE_DB.count()

	## Look at most popular reagents now
	zipsorted = sorted(counts.iteritems(), key = lambda x: x[1], reverse = True)
	sorted_ids = [x[0] for x in zipsorted]
	sorted_counts = [x[1] for x in zipsorted]
	ranks = range(1, len(sorted_counts) + 1)
	fig = plt.figure(figsize=(10,6), dpi = 300)
	ax = plt.gca()
	ax.scatter(ranks, sorted_counts, alpha = 0.5)
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.xlabel('Rank')
	plt.ylabel('Number of reaction instances')
	ax.axis([1, np.power(10, 1 + np.ceil(np.log10(max(ranks)))), \
	     0, max(sorted_counts) + 0.5])
	plt.title('Most popular solvents in Reaxys ({} instances)'.format(N_rxdids))
	plt.grid(True)
	N_top = 20
	for (i, _id) in enumerate(sorted_ids[:N_top]):
		label = names[_id]
		x = ranks[i]
		y = sorted_counts[i]

		plt.annotate(
			label, 
			xy = (x, y), xytext = (1, -(i - N_top + 0.5) / N_top),
			textcoords = 'axes fraction', ha = 'right', va = 'center',
			bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.8),
			arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
	fig.savefig(os.path.join(out_folder, 'solvents.png'))
	plt.show()

	with open(os.path.join(out_folder, 'solvents.csv'), 'w') as fid:
		fid.write('rank\t_id\tcount\tname\tsmiles\n')
		for (i, _id) in enumerate(sorted_ids):
			fid.write('{}\t{}\t{}\t{}\t{}\n'.format(
				i + 1, _id, counts[_id], names[_id], smiles[_id]
			))
	