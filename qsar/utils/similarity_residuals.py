'''
Look at the prediction residual against similarity scores
'''

import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

fname = '/home/ccoley/ML Chemistry/Make-It/makeit/models/'
fname += sys.argv[1]
mols = []
resids = []
with open(fname, 'r') as fid:
	for line in fid:
		smiles, resid = line.strip().split('\t')
		mol = Chem.MolFromSmiles(smiles)
		if not mol: 
			print(line)
		mols.append(mol)
		resids.append(float(resid))
[Chem.SanitizeMol(x) for x in mols]

#fps = [FingerprintMols.FingerprintMol(x) for x in mols]
fps = [AllChem.GetMorganFingerprint(x, 3) for x in mols]

N = len(mols)
similarities = np.zeros((N, N))
for i in range(len(mols)):
	print('Getting similarity score for {}'.format(i))
	for j in range(i + 1, len(mols)):
		#similarity = DataStructs.FingerprintSimilarity(fps[i], fps[j])
		similarity = DataStructs.DiceSimilarity(fps[i], fps[j])
		similarities[i, j] = similarity
		similarities[j, i] = similarity

		if similarity == 1:
			print(i)
			print(j)
			print(resids[i])
			print(resids[j])
			print(Chem.MolToSmiles(mols[i]))
			print(Chem.MolToSmiles(mols[j]))

sqresids = np.array([x ** 2 for x in resids])
absresids = np.array([abs(x) for x in resids])

meas = np.max(similarities, axis = 0)
plt.scatter(meas, absresids, alpha = 0.5)
plt.xlabel('Similarity score of most similar molecule')
plt.ylabel('Absolute error')
plt.title('Absolute error as a function of maximum similarity' + 
	      '\n{}'.format(fname.split('/')[-1].split(' ')[0]))
plt.grid(True)
plt.axis([0, max(meas), 0, max(absresids)])	
folder = os.path.dirname(fname)
plt.savefig(os.path.join(folder, 'max_similarity.png'), bbox_inches = 'tight')
plt.show()
plt.clf()

meas = np.sum(similarities, axis = 0) / (N - 1.)
plt.scatter(meas, absresids, alpha = 0.5)
plt.xlabel('Average similarity score in dataset')
plt.ylabel('Absolute error')
plt.title('Absolute error as a function of average similarity' + 
	      '\n{}'.format(fname.split('/')[-1].split(' ')[0]))
plt.grid(True)
plt.axis([0, max(meas), 0, max(absresids)])	
folder = os.path.dirname(fname)
plt.savefig(os.path.join(folder, 'avg_similarity.png'), bbox_inches = 'tight')
plt.show()
plt.clf()


similarities_copy = similarities.copy()
similarities_copy.sort()
meas = np.sum(similarities_copy[:, -5:], axis = 1) / 5.
print(meas.shape)
print(similarities_copy).shape


plt.scatter(meas, absresids, alpha = 0.5)
plt.xlabel('Average of top 5 similarity scores')
plt.ylabel('Absolute error')
plt.title('Absolute error as a function of top 5 similarities' + 
	      '\n{}'.format(fname.split('/')[-1].split(' ')[0]))
plt.grid(True)
plt.axis([0, max(meas), 0, max(absresids)])	
folder = os.path.dirname(fname)
plt.savefig(os.path.join(folder, 'avgmax5_similarity.png'), bbox_inches = 'tight')
plt.show()
plt.clf()

# Now take the sliding average from the previous graph
division = 0.1
xs = []; ys = []
for i in range(int(np.ceil(1/division))):
	x = i * division
	xs.append(x + division / 2)
	ys.append(absresids[(meas < x + division) * (meas > x)])
plt.scatter(xs, [np.mean(y) for y in ys], alpha = 1)
plt.errorbar(xs, [np.mean(y) for y in ys], [np.std(y) for y in ys], linestyle='None')
plt.xlabel('Binned - Average of top 5 similarity scores')
plt.ylabel('Averaged - Absolute error')
plt.title('Absolute error as a function of top 5 similarities' + 
	      '\n{}'.format(fname.split('/')[-1].split(' ')[0]))
plt.grid(True)	
folder = os.path.dirname(fname)
plt.savefig(os.path.join(folder, 'avgmax5_boxed_similarity.png'), bbox_inches = 'tight')
plt.show()
plt.clf()

means = [np.mean(y) for y in ys]
stds = [np.std(y) for y in ys]
with open(fname + '.data', 'w') as fid:
	fid.write('Center of bin (average of top 5 similarity scores)\tMAE\tAE(stdev)')
	for i in range(len(xs)):
		fid.write('{}\t{}\t{}\n'.format(xs[i], means[i], stds[i]))