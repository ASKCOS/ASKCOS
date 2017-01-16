# Import relevant packages
from __future__ import print_function
import time
from global_config import USE_STEREOCHEMISTRY
import numpy as np
import os
import sys
import argparse
import rdkit.Chem as Chem
import cPickle as pickle
import matplotlib
matplotlib.use('Agg')

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument('--data_tag', type = str, default = 'makeit/predict/data_edits_reaxys/reaxys',
		                help = 'Data file tag, default makeit/predict/data_edits_reaxys/reaxys')
	parser.add_argument('--lr', type = float, default = 0.01, 
						help = 'Learning rate, default 0.01')
	args = parser.parse_args()

	# Labels
	FROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'output')
	
	DATA_FPATH = '{}_data.pickle'.format(args.data_tag)


	from collections import defaultdict
	sizes_all = {
		'N_c':  defaultdict(int),
		'N_e1': defaultdict(int),
		'N_e2': defaultdict(int),
		'N_e3': defaultdict(int),
		'N_e4': defaultdict(int),
	}
	sizes_true = {
		'N_e1': defaultdict(int),
		'N_e2': defaultdict(int),
		'N_e3': defaultdict(int),
		'N_e4': defaultdict(int),
	}

	# Keep returning forever and ever
	k = 0
	with open(DATA_FPATH, 'rb') as fid:
			
		legend_data = pickle.load(fid) # first doc is legend

		# Pre-load indeces
		CANDIDATE_EDITS_COMPACT = legend_data['candidate_edits_compact']

				
		while True:
			try:
				doc = pickle.load(fid)
				k += 1
				if k % 1000 == 0: print('{} done'.format(k))
			except EOFError:
				print('End of file!')
				break

			sizes_all['N_c'][len(doc[CANDIDATE_EDITS_COMPACT])] += 1

			for (c, edit_string) in enumerate(doc[CANDIDATE_EDITS_COMPACT]):

				edit_string_split = edit_string.split(';')
				N_e1 = edit_string_split[0].count(',') + 1
				N_e2 = edit_string_split[1].count(',') + 1
				N_e3 = edit_string_split[2].count(',') + 1
				N_e4 = edit_string_split[3].count(',') + 1

				if c == 0:
					sizes_true['N_e1'][N_e1] += 1
					sizes_true['N_e2'][N_e2] += 1
					sizes_true['N_e3'][N_e3] += 1
					sizes_true['N_e4'][N_e4] += 1
				sizes_all['N_e1'][N_e1] += 1
				sizes_all['N_e2'][N_e2] += 1
				sizes_all['N_e3'][N_e3] += 1
				sizes_all['N_e4'][N_e4] += 1

	# print(sizes_all)
	# print(sizes_true)
	print('After looking at {} documents'.format(k))

	import matplotlib.pyplot as plt

	for i, tag in enumerate(['N_e1', 'N_e2', 'N_e3', 'N_e4']):
		plt.subplot(2, 2, i + 1)
		xs = range(20)
		ys = np.array([sizes_all[tag][x] for x in xs], dtype = float)
		ys /= np.sum(ys)
		ys_true = np.array([sizes_true[tag][x] for x in xs], dtype = float)
		ys_true /= np.sum(ys_true)
		plt.bar(xs, ys, 1.0, color='blue', alpha = 0.5, label = 'All candidates')
		plt.bar(xs, ys_true, 1.0, color='green', alpha = 0.5, label = 'True outcomes')
		plt.xlabel(tag)
		plt.ylabel('Normalized frequency')
		plt.legend()
		plt.yscale('log')
	plt.tight_layout()
	plt.savefig(os.path.join(os.path.dirname(DATA_FPATH), 'edit_analysis.png'))

	plt.clf()
	xs = range(0, 5000, 50)
	ys = np.zeros_like(np.array(xs), dtype = float)
	for (key, val) in sizes_all['N_c'].iteritems():
		for i, x in enumerate(xs):
			if key < x:
				ys[max(i-1, 0)] += val
				break
	ys /= np.sum(ys)
	plt.bar(xs, ys, 1.0, color='blue', alpha = 1.0)
	plt.xlabel('Number of candidates')
	plt.ylabel('Normalized frequency')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig(os.path.join(os.path.dirname(DATA_FPATH), 'edit_analysis_Nc.png'))