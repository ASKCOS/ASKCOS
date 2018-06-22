from makeit.prioritization.prioritizer import Prioritizer
import makeit.global_config as gc
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from makeit.utilities.buyable.pricer import Pricer
from makeit.utilities.io.logging import MyLogger
from keras.models import Sequential
from keras.layers import Dense, Lambda
from keras import backend as K
import math
import sys
import random
import os
import time
import os
import six; from six.moves import cPickle as pickle
from numpy import inf
mincost_prioritizer_loc = 'mincostprioritizer'


class MinCostPrecursorPrioritizer(Prioritizer):
    '''
    This is a standalone, importable MinCost model. Uses loaded keras model.
    '''

    def __init__(self, score_scale=10.0):
        
        self.vars = []
        self.FP_rad = 3
        self.score_scale = score_scale
        
        self._restored = False
        self.pricer = None
        self._loaded = False
 
    def load_model(self, FP_len=1024, input_layer=6144, hidden_layer=512, modelpath=""):
        
        self.FP_len = FP_len

        def last_layer(x):
            return (lambda x: 11 - 11 * K.exp(- K.abs(x / 1500000)))(x)

        model = Sequential()        
        model.add(Dense(input_layer, activation="relu", batch_input_shape=(None, self.FP_len) ))
        model.add(Dense(hidden_layer, activation="relu"))
        model.add(Dense(hidden_layer, activation="relu"))
        model.add(Dense(hidden_layer, activation="relu"))
        model.add(Dense(hidden_layer, activation="relu"))    
        model.add(Dense(1))
        model.add(Lambda(last_layer, output_shape = (None, 1)))
        if modelpath == "":
            modelpath = gc.MinCost_Prioritiaztion['trained_model_path']
        weights = model.load_weights(modelpath)
        self.model = model 
     
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

    def get_price(self, smi):
        ppg = self.pricer.lookup_smiles(smi, alreadyCanonical=True)
        if ppg:
            return 0.0
        else:
            return None 

    def get_priority(self, retroProduct, **kwargs):
        mode = kwargs.pop('mode', gc.max)
        if not self._loaded:
            self.load_model()

        if not isinstance(retroProduct, str):
            scores = []
            for smiles in retroProduct.smiles_list:
                scores.append(self.get_score_from_smiles(smiles))
            return sum(scores)
        else:
            return self.get_score_from_smiles(retroProduct)
        if not retroProduct:
            return inf

    def get_score_from_smiles(self, smiles):
        # Check buyable
        ppg = self.pricer.lookup_smiles(smiles, alreadyCanonical=True)
        if ppg:
            return 0.0 #ppg / 100.
        
        fp = np.array((self.smi_to_fp(smiles)), dtype=np.float32)
        if sum(fp) == 0:
            cur_score = 0.
        else:
            cur_score = self.model.predict(fp.reshape((1,self.FP_len)))[0][0]
        return cur_score

if __name__ == '__main__':
    
    model = MinCostPrecursorPrioritizer()
    model.load_model()
    smis = ['CC(=O)N1C=C(C=C2N=C(N(N=CC3C=CC=CC=3)C(C)=O)N(C(C)=O)C2=O)C2=CC=CC=C21', 'CCCNc1ccccc1']
    for smi in smis:
        sco = model.get_priority(smi)
        print('{} <--- {}'.format(sco, smi))
