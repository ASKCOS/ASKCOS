'''
This is a standalone, importable SCScorer model. It does not have tensorflow as a
dependency and is a more attractive option for deployment. The calculations are 
fast enough that there is no real reason to use GPUs (via tf) instead of CPUs (via np)
'''

import math, sys, random, os
import numpy as np
import time
import rdkit.Chem as Chem 
import rdkit.Chem.AllChem as AllChem

import os 
project_root = os.path.dirname(os.path.dirname(__file__))

score_scale = 5.0
min_separation = 0.25

FP_len = 1024
FP_rad = 2

def mol_to_fp(mol, radius=FP_rad, nBits=FP_len):
    if mol is None:
        return np.zeros((nBits,), dtype=np.float32)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits, 
        useChirality=True), dtype=np.bool)

def smi_to_fp(smi, radius=FP_rad, nBits=FP_len):
    if not smi:
        return np.zeros((nBits,), dtype=np.float32)
    return mol_to_fp(Chem.MolFromSmiles(smi), radius, nBits)

def sigmoid(x):
  return 1 / (1 + math.exp(-x))

def softmax(x):
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum()

class RetroTempPrioritizer():
    def __init__(self, score_scale=score_scale):
        self.vars = []
        self.score_scale = score_scale
        self._restored = False

    def restore(self, weight_path=os.path.join(project_root, 'models', '5d4M_Reaxys', 'model.ckpt-105660.as_numpy.pickle')):
        import cPickle as pickle
        with open(weight_path, 'rb') as fid:
            self.vars = pickle.load(fid)
        print('Restored variables from {}'.format(weight_path))
        self._restored = True
        return self

    def apply(self, x):
        if not self._restored:
            raise ValueError('Must restore model weights!')
        # Each pair of vars is a weight and bias term
        for i in range(0, len(self.vars), 2):
            last_layer = (i == len(self.vars)-2)
            W = self.vars[i] 
            b = self.vars[i+1]
            x = np.matmul(x, W) + b
            if not last_layer:
                x = x * (x > 0) # ReLU
        return x

    def get_topk_from_smi(self, smi='', k=100):
        if not smi:
            return []
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            return []
        return self.get_topk_from_mol(mol, k=k)
        
    def get_topk_from_mol(self, mol, k=100):
        fp = mol_to_fp(mol).astype(np.float32)
        cur_scores = self.apply(fp)
        indices = list(cur_scores.argsort()[-k:][::-1])
        cur_scores.sort()
        probs = softmax(cur_scores)
        return probs[-k:][::-1], indices


if __name__ == '__main__':
    model = RetroTempPrioritizer()    
    model.restore(os.path.join(project_root, 'models', '5d4M_Reaxys', 'model.ckpt-105660.as_numpy.pickle'))
    smis = ['CCCOCCC', 'CCCNc1ccccc1']
    for smi in smis:
        probs, lst = model.get_topk_from_smi(smi)
        print('{} -> {}'.format(smi, lst))