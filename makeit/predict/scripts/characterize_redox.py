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
	REDOX_DB = db['redox']
	CHEMICAL_DB = db['chemicals']
	
	counts = []
	powers = []
	labels = []
	ids = []
	refs = set()
	N_rxdids = 0

	with open(os.path.join(out_folder, 'redox.dat'), 'w') as fid:
		for doc in REDOX_DB.find({}, ['_id', 'SMILES', 'count', 'ox_power', 'references']):
			fid.write('{}\t{}\t{}\t{}\n'.format(doc['_id'], doc['SMILES'], doc['count'], doc['ox_power']))

			counts.append(doc['count'])
			powers.append(doc['ox_power'])
			labels.append(doc['SMILES'])
			ids.append(doc['_id'])
			refs |= set(doc['references'])
			
	N_rxdids = len(refs)

	weighted_powers = [np.log10(count) * power for (count, power) in zip(counts, powers)]

	# Oxidizing agent scatterplot
	fig = plt.figure(figsize=(6,4), dpi = 300)
	ax = plt.gca()
	ax.scatter(powers, counts, alpha = 0.5)
	ax.set_yscale('log')
	plt.xlabel('Average oxidizing power')
	plt.ylabel('Count (number of examples)')
	plt.title('Redox power of reagents in Reaxys ({} examples)'.format(N_rxdids))
	plt.grid(True)
	ax.axis([min(powers) - 1, max(powers) + 1, \
	     0.8, np.power(10, 0.2 + np.ceil(np.log10(max(counts))))])
	fig.savefig(os.path.join(out_folder, 'redox.png'))

	## Now look at oxidizing specificially
	zipsorted = sorted(zip(weighted_powers, labels, ids), key = lambda x: x[0])
	weighted_powers = np.array([x[0] for x in zipsorted])
	labels = [x[1] for x in zipsorted]
	ids = [x[2] for x in zipsorted]
	ranks = range(1, len(weighted_powers) + 1)
	ranks.reverse()
	fig = plt.figure(figsize=(8,4), dpi = 300)
	ax = plt.gca()
	ax.scatter(ranks, weighted_powers, alpha = 0.5)
	ax.set_xscale('log')
	plt.xlabel('Rank')
	plt.ylabel('Average power * log10(count)')
	ax.axis([1, np.power(10, 1 + np.ceil(np.log10(max(ranks)))), \
	     0, max(weighted_powers) + 0.5])
	plt.title('Perceived oxidizing agents in Reaxys ({} examples)'.format(N_rxdids))
	plt.grid(True)
	N_top = 8
	for (i, (label, x, y, _id)) in enumerate(zip(labels, ranks, weighted_powers, ids)[-N_top:]):

		chem_doc = CHEMICAL_DB.find_one({'_id': _id}, ['IDE_CN'])
		label = str(chem_doc['IDE_CN'])

		plt.annotate(
			label, 
			xy = (x, y), xytext = (1, 1 + (i - N_top + 0.5) / N_top),
			textcoords = 'axes fraction', ha = 'right', va = 'center',
			bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.8),
			arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
	fig.savefig(os.path.join(out_folder, 'redox - ox.png'))

	## Now look at oxidizing specificially
	zipsorted = sorted(zip(weighted_powers, labels, ids), key = lambda x: -x[0])
	weighted_powers = np.array([x[0] for x in zipsorted])
	labels = [x[1] for x in zipsorted]
	ids = [x[2] for x in zipsorted]
	ranks = range(1, len(weighted_powers) + 1)
	ranks.reverse()
	fig = plt.figure(figsize=(8,4), dpi = 300)
	ax = plt.gca()
	ax.scatter(ranks, weighted_powers, alpha = 0.5)
	ax.set_xscale('log')
	plt.xlabel('Rank')
	plt.ylabel('Average power * log10(count)')
	ax.axis([1, np.power(10, 1 + np.ceil(np.log10(max(ranks)))), \
	     min(weighted_powers) - 0.5, 0])
	plt.title('Perceived reducing agents in Reaxys ({} examples)'.format(N_rxdids))
	plt.grid(True)
	N_top = 8
	for (i, (label, x, y, _id)) in enumerate(zip(labels, ranks, weighted_powers, ids)[-N_top:][::-1]):

		chem_doc = CHEMICAL_DB.find_one({'_id': _id}, ['IDE_CN'])
		label = str(chem_doc['IDE_CN'])

		plt.annotate(
			label, 
			xy = (x, y), xytext = (1, (i - N_top + 0.5) / N_top),
			textcoords = 'axes fraction', ha = 'right', va = 'center',
			bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.8),
			arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
	fig.savefig(os.path.join(out_folder, 'redox - red.png'))