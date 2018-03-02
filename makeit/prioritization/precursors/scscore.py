from makeit.prioritization.prioritizer import Prioritizer
import makeit.global_config as gc
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from makeit.utilities.buyable.pricer import Pricer
from makeit.utilities.io.logging import MyLogger
import math
import sys
import random
import os
import time
import os
import cPickle as pickle
from numpy import inf
scscore_prioritizer_loc = 'scscoreprioritizer'


class SCScorePrecursorPrioritizer(Prioritizer):
    '''
    This is a standalone, importable SCScorecorer model. It does not have tensorflow as a
    dependency and is a more attractive option for deployment. The calculations are 
    fast enough that there is no real reason to use GPUs (via tf) instead of CPUs (via np)
    '''

    def __init__(self, score_scale=5.0):
        self.vars = []
        self.FP_rad = 2
        self.score_scale = score_scale
        self._restored = False
        self.pricer = None
        self._loaded = False

    def load_model(self, FP_len=1024, model_tag='1024bool'):
        self.FP_len = FP_len
        if model_tag != '1024bool' and model_tag != '1024uint8' and model_tag != '2048bool':
            MyLogger.print_and_log(
                'Non-existent SCScore model requested: {}. Using "1024bool" model'.format(model_tag), scscore_prioritizer_loc, level=2)
            model_tag = '1024bool'
        filename = 'trained_model_path_'+model_tag
        with open(gc.SCScore_Prioritiaztion[filename], 'rb') as fid:
            self.vars = pickle.load(fid)
        if gc.DEBUG:
            MyLogger.print_and_log('Loaded synthetic complexity score prioritization model from {}'.format(
            gc.SCScore_Prioritiaztion[filename]), scscore_prioritizer_loc)

        if 'uint8' in gc.SCScore_Prioritiaztion[filename]:
            def mol_to_fp(mol):
                if mol is None:
                    return np.array((self.FP_len,), dtype=np.uint8)
                fp = AllChem.GetMorganFingerprint(
                    mol, self.FP_rad, useChirality=True)  # uitnsparsevect
                fp_folded = np.zeros((self.FP_len,), dtype=np.uint8)
                for k, v in fp.GetNonzeroElements().iteritems():
                    fp_folded[k % self.FP_len] += v
                return np.array(fp_folded)
        else:
            def mol_to_fp(mol):
                if mol is None:
                    return np.zeros((self.FP_len,), dtype=np.float32)
                return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, self.FP_rad, nBits=self.FP_len,
                                                                      useChirality=True), dtype=np.bool)
        self.mol_to_fp = mol_to_fp

        self.pricer = Pricer()
        self.pricer.load_from_file()
        self._restored = True
        self._loaded = True

    def smi_to_fp(self, smi):
        if not smi:
            return np.zeros((self.FP_len,), dtype=np.float32)
        return self.mol_to_fp(Chem.MolFromSmiles(smi))

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
                x = x * (x > 0)  # ReLU
        x = 1 + (self.score_scale - 1) * sigmoid(x)
        return x

    def get_priority(self, retroProduct, **kwargs):
        mode = kwargs.get('mode', gc.max)
        if not self._loaded:
            self.load_model()

        if not isinstance(retroProduct, str):
            scores = []
            for smiles in retroProduct.smiles_list:
                scores.append(self.get_score_from_smiles(smiles))
            return -self.merge_scores(scores, mode = mode)
        else:
            return -self.get_score_from_smiles(retroProduct)
        if not retroProduct:
            return -inf

    def merge_scores(self, list_of_scores, mode=gc.max):
        if mode == gc.mean:
            return np.mean(list_of_scores)
        elif mode == gc.geometric:
            return np.power(np.prod(list_of_scores), 1.0/len(list_of_scores))
        elif mode == gc.pow8:
            pow8 = []
            for score in list_of_scores:
                pow8.append(8**score)
            return np.sum(pow8)
        else:
            return np.max(list_of_scores)

    def get_score_from_smiles(self, smiles):
        # Check buyable
        ppg = self.pricer.lookup_smiles(smiles, alreadyCanonical=True)
        if ppg:
            return ppg / 100.
        
        fp = np.array((self.smi_to_fp(smiles)), dtype=np.float32)
        if sum(fp) == 0:
            cur_score = 0.
        else:
            # Run
            cur_score = self.apply(fp)
        return cur_score


def sigmoid(x):
    return 1 / (1 + math.exp(-x))

if __name__ == '__main__':
    model = SCScorer()
    model.load_model(model_tag='1024bool')
    smis = ['CCCOCCC', 'CCCNc1ccccc1']
    for smi in smis:
        sco = model.get_priority(smi)
        print('{} <--- {}'.format(sco, smi))

    model = SCScorer()
    model.load_model(model_tag='2048bool', FP_len=2048)
    smis = ['CCCOCCC', 'CCCNc1ccccc1']
    for smi in smis:
        sco = model.get_priority(smi)
        print('{} <--- {}'.format(sco, smi))

    model = SCScorer()
    model.load_model(model_tag='1024uint8')
    smis = ['CCCOCCC', 'CCCNc1ccccc1']
    for smi in smis:
        sco = model.get_priority(smi)
        print('{} <--- {}'.format(sco, smi))
