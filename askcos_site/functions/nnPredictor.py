import os
import argparse
import time
import cPickle as pickle
import numpy as np
import random
import matplotlib.pyplot as plt
from termcolor import colored
from collections import Counter, defaultdict
from sklearn.neighbors import NearestNeighbors as NN
from sklearn.externals import joblib
# from draw import ReactionStringToImage
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

# DATABASE
from pymongo import MongoClient  # mongodb plugin

client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaxys']
REACTION_DB = db['reactions']
INSTANCE_DB = db['instances']
CHEMICAL_DB = db['chemicals']
SOLVENT_DB = db['solvents']


# # Load all the instance IDs from the test model
# project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# model_dir = os.path.join(project_dir, 'testModel')
# rxd_ids = []
# rxn_ids = []
# with open(os.path.join(model_dir, '1650000-1699999_10NN_20000SRR_info.txt'), 'r') as infile:
#     rxn_ids.append(infile.readlines()[1:])  # a list of str(rxn_ids) with '\n'
# for id in rxn_ids[0]:
#     rxd_ids.append(id.replace('\n', ''))
# # Load the test NN model
# lshf_nn = joblib.load(os.path.join(model_dir, '1650000-1699999_10NN_20000SRR_lshf.pickle'))


# # Load the full model
# model_dir = '/home/yrwang/askcos/2MRxnModel'
# rxd_ids = []
# rxn_ids = []
# with open(os.path.join(model_dir, 'fpNN-10_2MRxn_info.txt'), 'r') as infile:
#     rxn_ids.append(infile.readlines()[1:])  # a list of str(rxn_ids) with '\n'
# for id in rxn_ids[0]:
#     rxd_ids.append(id.replace('\n', ''))
# # Load the NN model
# lshf_nn = joblib.load(os.path.join(model_dir, 'fpNN-10_2MRxn_lshf.pickle'))


def string_or_range_to_float(text):
    try:
        return float(text)
    except Exception as e:
        if text.count('-') == 1:  # 20 - 30
            try:
                x = text.split('-')
                return (float(x[0]) + float(x[1])) / 2.0
            except Exception as e:
                print(e)
        elif text.count('-') == 2:  # -20 - 0
            try:
                x = text.split('-')
                return (-float(x[0]) + float(x[1])) / 2.0
            except Exception as e:
                print(e)
        elif text.count('-') == 3:  # -20 - -10
            try:
                x = text.split('-')
                return (-float(x[0]) - float(x[1])) / 2.0
            except Exception as e:
                print(e)
        else:
            print(e)
    return None


def instance_rxn_condition(INSTANCE_DB, int_id):
    """Return the reaction conditions of a particular instance with known instance ID"""
    doc = INSTANCE_DB.find_one({'_id': int_id})

    context_info = ''
    # Temp
    T = string_or_range_to_float(doc['RXD_T'])  # It can be None

    # Solvent(s)
    solvent = ''
    context_info += 'solv:'
    for xrn in doc['RXD_SOLXRN']:
        slvt = CHEMICAL_DB.find_one({'_id': xrn})
        if not slvt:
            print('########## COULD NOT FIND SOLVENT {} ###########'.format(xrn))
            continue
        smi_or_name = str(slvt['SMILES'])
        if not smi_or_name:
            smi_or_name = str(slvt['IDE_CN'])
        solvent += smi_or_name + '.'
        context_info += str(slvt['IDE_CN']) + '(' + str(slvt['SMILES']) + ')' + ','

    # Reagents
    reagent = ''
    context_info += 'rgt:'
    for xrn in doc['RXD_RGTXRN']:
        rgt = CHEMICAL_DB.find_one({'_id': xrn})
        if not rgt:
            print('########## COULD NOT FIND REAGENT {} ###########'.format(xrn))
            continue
        smi_or_name = str(rgt['SMILES'])
        if not smi_or_name:
            smi_or_name = str(rgt['IDE_CN'])
        reagent += smi_or_name + '.'
        context_info += str(rgt['IDE_CN']) + '(' + str(rgt['SMILES']) + ')' + ','

    # Catalysts
    catalyst = ''
    context_info += 'cat:'
    for xrn in doc['RXD_CATXRN']:
        cat = CHEMICAL_DB.find_one({'_id': xrn})
        if not cat:
            print('########## COULD NOT FIND CATALYST {} ###########'.format(xrn))
            continue
        smi_or_name = str(cat['SMILES'])
        if not smi_or_name:
            smi_or_name = str(cat['IDE_CN'])
        catalyst += smi_or_name + '.'
        context_info += str(cat['IDE_CN']) + '(' + str(cat['SMILES']) + ')' + ','

    # Time, yield
    rxn_time = doc['RXD_TIM']
    rxn_yield = doc['RXD_NYD']

    context_info += 'T:{}C'.format(T) + ',t:{}min'.format(rxn_time) + ',y:{}%'.format(rxn_yield)

    return context_info, T, rxn_time, rxn_yield, solvent, reagent, catalyst


def create_rxn_Morgan2FP(rsmi, psmi, fpsize=1024, useFeatures=True):
    """Create a rxn Morgan (r=2) fingerprint as bit vector from SMILES string lists of reactants and products"""
    # Modified from Schneider's code (2014)

    rfp = None
    for react in rsmi:
        mol = Chem.MolFromSmiles(react)
        try:
            fp = np.array(
                AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits=fpsize, useFeatures=useFeatures))
        except Exception as e:
            print("Cannot build reactant fp due to {}".format(e))
        if rfp is None:
            rfp = fp
        else:
            rfp += fp
    pfp = None
    for product in psmi:
        mol = Chem.MolFromSmiles(product)
        try:
            fp = np.array(
                AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits=fpsize, useFeatures=useFeatures))
        except Exception as e:
            print("Cannot build product fp due to {}".format(e))
        if pfp is None:
            pfp = fp
        else:
            pfp += fp
    if pfp is not None and rfp is not None:
        pfp -= rfp
    return pfp


def rxn_condition_predictor_amongNN(INSTANCE_DB, dists, idx, rxd_ids, num_c=1, dist_limit=0.3, outputString=True):
    """Reaction condition recommendation among the top 10 NN, maximum 2 recommendations

        dists, idx: np.array output from NearestNeighbor model for one rxn (10L, )
        The second recommendation is based on the most popular reagents within dist_limit
    """
    if num_c > 2:
        print('No more than 2 recommendations!')
        return None
    else:
        int_id = rxd_ids[idx[0]]
        context1, T1, t1, y1, slvt1, rgt1, cat1 = instance_rxn_condition(INSTANCE_DB, int_id)
        if num_c == 1 or dists[0] > dist_limit:
            print('Conditions for Top1 NN is used')
            if dists[0] > dist_limit:
                print('No neighbor is found within a cosine distance of {}'.format(dist_limit))

            if outputString:
                return context1
            else:
                return T1, t1, y1, slvt1, rgt1, cat1
        else:
            n_id = []
            for i, dist in enumerate(dists[1:]):
                if dist <= dist_limit:
                    n_id.append(idx[i + 1])
            if len(n_id) < 3:
                int_id2 = rxd_ids[idx[1]]
                context2, T2, t2, y2, slvt2, rgt2, cat2 = instance_rxn_condition(INSTANCE_DB, int_id2)
                if len(n_id) == 0:
                    print('No second neighbor is found within a cosine distance of {}'.format(dist_limit))
            else:
                rt = []
                for ids in n_id:
                    rr = instance_rxn_condition(INSTANCE_DB, rxd_ids[ids])[5]
                    rt.append(rr)
                rts = tuple(rt)
                rgt_counter = Counter(rts)
                rgt = rgt_counter.most_common(1)[0][0]
                for r in range(len(rt)):
                    if rt[r] == rgt:
                        int_id2 = rxd_ids[n_id[r]]
                        context2, T2, t2, y2, slvt2, rgt2, cat2 = instance_rxn_condition(INSTANCE_DB, int_id2)
                        break

            if outputString:
                return context1, context2
            else:
                return [T1, t1, y1, slvt1, rgt1, cat1], [T2, t2, y2, slvt2, rgt2, cat2]


def n_rxn_condition(INSTANCE_DB, n, dists, idx, rxd_ids, dist_limit=0.3):
    """Reaction condition list from the top 10 NN

    :param n: int, the number of nearest neighbors to extract rxn conditions from, n <= 10 here
    :param dists, idx: np.array output from NearestNeighbor model for one rxn (10L, )
    :param rxd_ids: the list of instance IDs for all the instances in the database
    :return: A list of reaction conditions [(temp, time, yield, solvents, reagents, cats), (), ]
    """
    if n > int(idx.shape[0]):
        print('More rxn condition options requested than the number of NN, n is set to {}'.format(idx.shape[1]))
    if dists[0] > dist_limit:
        print('No neighbor is found within a cosine distance of {}'.format(dist_limit))

    contexts = []
    for i, rid in enumerate(idx):
        if i >= n:
            break
        contexts.append(instance_rxn_condition(INSTANCE_DB, rxd_ids[rid])[1:]) # T, t, y, slvt, rgt, cat
    return contexts


class NNConditionPredictor():
    """Reaction condition predictor based on Nearest Neighbor method"""

    def __init__(self, nn_model=None, rxn_ids=None, INSTANCE_DB=None):
        self.nnModel = nn_model
        self.rxn_ids = rxn_ids
        self.num_cond = 1
        self.dist_limit = 0.3
        self.outputString = True
        self.INSTANCE_DB = INSTANCE_DB

    def load_predictor(self, userInput):
        self.num_cond = userInput['num_cond']
        self.dist_limit = userInput['dist_limit']
        self.outputString = userInput['outputString']

    def step_condition(self, rxn):
        """Reaction condition recommendation for a reaction

            rxn: [list of reactant SMILES strings, list of product SMILES strings]
            return: if outputString = True, context info as a string a tuple of strings; else lists of condition lists
        """
        rxn_fp = create_rxn_Morgan2FP(rxn[0], rxn[1], fpsize=1024, useFeatures=True)
        dists, ids =self.nnModel.kneighbors(rxn_fp)  # (1L, 10L)
        return rxn_condition_predictor_amongNN(self.INSTANCE_DB, dists[0], ids[0], rxd_ids=self.rxn_ids, num_c=self.num_cond,
                                               dist_limit=self.dist_limit, outputString=self.outputString)

    def step_n_conditions(self, n, rxn):
        """n reaction condition recommendations for a reaction

            rxn: [list of reactant SMILES strings, list of product SMILES strings]
            return: lists of condition tuples
        """
        rxn_fp = create_rxn_Morgan2FP(rxn[0], rxn[1], fpsize=1024, useFeatures=True)
        dists, ids =self.nnModel.kneighbors(rxn_fp)  # (1L, 10L)
        return n_rxn_condition(INSTANCE_DB, n, dists[0], ids[0], rxd_ids=self.rxn_ids, dist_limit=self.dist_limit)

    def path_condition(self, path):
        """Reaction condition recommendation for a reaction path with multiple reactions

            path: [[reactants, products], [reactants, products], ...]
            return: if outputString = True, context info as a lists of strings or tuple of strings;
            else lists of lists of condition lists
        """
        rxn_fps = []
        for react, product in path:
            rxn_fps.append(create_rxn_Morgan2FP(react, product, fpsize=1024, useFeatures=True))
        rxn_fps = np.array(rxn_fps)
        dists, ids = self.nnModel.kneighbors(rxn_fps)  # (nL, 10L)
        contexts = []
        for i, dist in enumerate(dists):
            contexts.append(rxn_condition_predictor_amongNN(INSTANCE_DB, dist, ids[i], rxd_ids=self.rxn_ids, num_c=self.num_cond,
                                                            dist_limit=self.dist_limit, outputString=self.outputString))
        return contexts
