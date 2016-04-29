'''
Look at the prediction residual against similarity scores
'''

import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import matplotlib.pyplot as plt
import numpy as np

fname = '/home/ccoley/ML Chemistry/Make-It/makeit/models/'
fname += 'Ab-oct-1a/Ab-oct-1a all.test'
mols = []
resids = []
with open(fname, 'r') as fid:
	for line in fid:
		smiles, resid = line.strip().split('\t')
		mol = Chem.MolFromSmiles(smiles)
		if not mol: 
			print(line)
		mols.append(mol)
		resids.append(resid)
[Chem.SanitizeMol(x) for x in mols]

#fps = [FingerprintMols.FingerprintMol(x) for x in mols]
fps = [AllChem.GetMorganFingerprint(x, 5) for x in mols]

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

meas = [0] * N 
for i in range(N):
	# avg
	#meas = np.sum(similarities, axis = 0) / (N - 1.)

	# max
	meas = np.max(similarities, axis = 0)

plt.scatter(meas, resids, alpha = 0.5)
plt.xlabel('Tanimoto score of most similar molecule')
plt.ylabel('Residual (Ab-oct-1a)')
plt.title('Residual as a function of similarity')
plt.grid(True)
plt.axis([0, 1, -2, 2])	
#plt.savefig(test_fpath + ' {}.png'.format(set_label), bbox_inches = 'tight')
plt.show()
plt.clf()