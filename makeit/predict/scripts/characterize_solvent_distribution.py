# Import relevant packages
from makeit.utils.chemnet_connect import *      # mongodb connection, gets 'chemicals', 'reactions'
import numpy as np     	      	   # for simple calculations
import matplotlib.pyplot as plt    # for visualization
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import os                          # for saving
from tqdm import tqdm
import cPickle as pickle
from collections import defaultdict

	

if __name__ == '__main__':
	out_folder = os.path.join(os.path.dirname(__file__), 'output')

	# DATABASE
	from pymongo import MongoClient    # mongodb plugin
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaxys']
	INSTANCE_DB = db['instances']
	CHEMICAL_DB = db['chemicals']

	if os.path.isfile(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), 'solvents.pickle')):
		with open(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), 'solvents.pickle'), 'rb') as fid:
			solvents = pickle.load(fid)
	else:
		solvents = defaultdict(int) # dict of counts
		for doc in INSTANCE_DB.find({}, ['_id', 'RXD_SOLXRN']):
			for xrn in doc['RXD_SOLXRN']:
				solvents[xrn] += 1

	counts = solvents.values()
	solvent_ids = solvents.keys()
	smiles = [CHEMICAL_DB.find_one({'_id': _id})['SMILES'] if CHEMICAL_DB.find_one({'_id': _id}) else _id for _id in solvent_ids]
	total_count = sum(counts)

	zipsorted = sorted(zip(counts, smiles), key = lambda x: x[0])
	counts = np.array([x[0] for x in zipsorted])
	labels = [x[1] for x in zipsorted]
	probs = counts * 1.0 / np.sum(counts)
	ranks = range(1, len(probs) + 1)
	ranks.reverse()

	# Oxidizing agent scatterplot
	fig = plt.figure(figsize=(6,4), dpi = 300)
	ax = plt.gca()
	ax.scatter(ranks, probs, alpha = 0.5)
	ax.set_yscale('log')
	ax.set_xscale('log')
	plt.xlabel('Rank')
	plt.ylabel('Probability')
	ax.axis([1, np.power(10, np.ceil(np.log10(max(ranks)))), \
	     np.power(10, np.floor(np.log10(min(probs)))), \
	     np.power(10, np.ceil(np.log10(max(probs))))])
	plt.title('{} solvent records in {} Reaxys instances'.format(total_count, INSTANCE_DB.count()))
	plt.grid(True)
	N_top = 10
	for (i, (label, x, y)) in enumerate(zip(labels, ranks, probs)[-N_top:]):
		plt.annotate(
			label, 
			xy = (x, y), xytext = (1, 1 + (i - N_top + 0.5) / N_top),
			textcoords = 'axes fraction', ha = 'right', va = 'center',
			bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.8),
			arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
	fig.savefig(os.path.join(out_folder, 'solvents.png'))



	# Now label with name

	counts = solvents.values()
	solvent_ids = solvents.keys()
	names = [CHEMICAL_DB.find_one({'_id': _id})['IDE_CN'] if CHEMICAL_DB.find_one({'_id': _id}) else _id for _id in solvent_ids]
	total_count = sum(counts)

	zipsorted = sorted(zip(counts, names), key = lambda x: x[0])
	counts = np.array([x[0] for x in zipsorted])
	labels = [x[1] for x in zipsorted]
	probs = counts * 1.0 / np.sum(counts)
	ranks = range(1, len(probs) + 1)
	ranks.reverse()

	fig = plt.figure(figsize=(6,4), dpi = 300)
	ax = plt.gca()
	ax.scatter(ranks, probs, alpha = 0.5)
	ax.set_yscale('log')
	ax.set_xscale('log')
	plt.xlabel('Rank')
	plt.ylabel('Probability')
	ax.axis([1, np.power(10, np.ceil(np.log10(max(ranks)))), \
	     np.power(10, np.floor(np.log10(min(probs)))), \
	     np.power(10, np.ceil(np.log10(max(probs))))])
	plt.title('{} solvent records in {} Reaxys instances'.format(total_count, INSTANCE_DB.count()))
	plt.grid(True)
	N_top = 10
	for (i, (label, x, y)) in enumerate(zip(labels, ranks, probs)[-N_top:]):
		plt.annotate(
			label, 
			xy = (x, y), xytext = (1, 1 + (i - N_top + 0.5) / N_top),
			textcoords = 'axes fraction', ha = 'right', va = 'center',
			bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.8),
			arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
	fig.savefig(os.path.join(out_folder, 'solvents_name.png'))

	with open(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), 'solvents.pickle'), 'wb') as fid:
		pickle.dump((solvents), fid, pickle.HIGHEST_PROTOCOL)
	with open(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), 'solvents.csv'), 'w') as fid:
		for (label, count) in zip(labels, counts):
			fid.write('{}\t{}\n'.format(label, count))