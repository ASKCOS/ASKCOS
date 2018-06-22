from makeit.utilities.fastfilter_utilities import Highway_self, pos_ct, true_pos, real_pos, set_keras_backend
from makeit.utilities.fingerprinting import create_rxn_Morgan2FP_separately
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from makeit.interfaces.scorer import Scorer
import numpy as np
import csv
from pymongo import MongoClient
from tqdm import tqdm
from keras.models import load_model
from keras import backend as K
import makeit.global_config as gc
from makeit.utilities.io.logging import MyLogger
import os
fast_filter_loc = 'fast_filter'


class FastFilterScorer(Scorer):
    def __init__(self):
        self.model = None

    def set_keras_backend(self, backend):
        if K.backend() != backend:
            os.environ['KERAS_BACKEND'] = backend
            reload(K)
            assert K.backend() == backend

    def load(self, model_path):
        MyLogger.print_and_log('Starting to load fast filter', fast_filter_loc)
        self.model = load_model(model_path, custom_objects={
                                'Highway_self': Highway_self, 'pos_ct': pos_ct, 'true_pos': true_pos, 'real_pos': real_pos})
        self.model._make_predict_function()
        MyLogger.print_and_log('Done loading fast filter', fast_filter_loc)

    def evaluate(self, reactant_smiles, target, **kwargs):
        # Strip chirality
        # rmol = Chem.MolFromSmiles(reactant_smiles)
        # pmol = Chem.MolFromSmiles(target)
        # reactant_smiles = Chem.MolToSmiles(rmol, False)
        # target = Chem.MolToSmiles(pmol, False)

        [pfp, rfp] = create_rxn_Morgan2FP_separately(
            reactant_smiles, target, rxnfpsize=2048, pfpsize=2048, useFeatures=False)
        pfp = np.asarray(pfp, dtype='float32')
        rfp = np.asarray(rfp, dtype='float32')
        rxnfp = pfp - rfp

        score = self.model.predict(
            [pfp.reshape(1, 2048), rxnfp.reshape(1, 2048)])
        outcome = {'smiles': target,
                   'template_ids': [],
                   'num_examples': 0
                   }
        all_outcomes = []
        all_outcomes.append([{'rank': 1.0,
                              'outcome': outcome,
                              'score': float(score[0][0]),
                              'prob': float(score[0][0]),
                              }])
        return all_outcomes

    def filter_with_threshold(self, reactant_smiles, target, threshold):
        [pfp, rfp] = create_rxn_Morgan2FP_separately(
            reactant_smiles, target, rxnfpsize=2048, pfpsize=2048, useFeatures=False)
        pfp = np.asarray(pfp, dtype='float32')
        rfp = np.asarray(rfp, dtype='float32')
        rxnfp = pfp - rfp

        score = self.model.predict([pfp.reshape(1, 2048), rxnfp.reshape(1, 2048)])
        filter_flag = (score > threshold)
        return filter_flag, float(score)


if __name__ == "__main__":

    ff = FastFilterScorer()
    ff.load(model_path=gc.FAST_FILTER_MODEL['trained_model_path'])
    score = ff.evaluate('CCO.CC(=O)O', 'CCOC(=O)C')
    print(score)
    score = ff.evaluate('[CH3:1][C:2](=[O:3])[O:4][CH:5]1[CH:6]([O:7][C:8]([CH3:9])=[O:10])[CH:11]([CH2:12][O:13][C:14]([CH3:15])=[O:16])[O:17][CH:18]([O:19][CH2:20][CH2:21][CH2:22][CH2:23][CH2:24][CH2:25][CH2:26][CH2:27][CH2:28][CH3:29])[CH:30]1[O:31][C:32]([CH3:33])=[O:34].[CH3:35][O-:36].[CH3:38][OH:39].[Na+:37]', 'CCCCCCCCCCOC1OC(CO)C(O)C(O)C1O')
    print(score)
    score = ff.evaluate(
        'CNC.Cc1ccc(S(=O)(=O)OCCOC(c2ccccc2)c2ccccc2)cc1', 'CN(C)CCOC(c1ccccc1)c2ccccc2')
    print(score)

    flag, sco = ff.filter_with_threshold('CCO.CC(=O)O', 'CCOC(=O)C', 0.75)
    print(flag)
    print(sco)
