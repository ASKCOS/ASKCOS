import os
import cPickle as pickle
import numpy as np
import random
# from sklearn.neighbors import NearestNeighbors as NN
from sklearn.externals import joblib
from rdkit import Chem
from rdkit.Chem import AllChem

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


def create_rxn_Morgan2FP(rxn_smiles, fpsize=1024, useFeatures=True):
    """Create a rxn Morgan (r=2) fingerprint as bit vector from a reaction SMILES string

        Modified from Schneider's code (2014)"""

    rsmi = rxn_smiles.split('>')[0].split('.')
    psmi = rxn_smiles.split('>')[2].split('.')

    rfp = None
    pfp = None
    for react in rsmi:
        mol = Chem.MolFromSmiles(react)
        fp = np.array(
            AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits=fpsize, useFeatures=useFeatures))
        if rfp is None:
            rfp = fp
        else:
            rfp += fp
    for product in psmi:
        mol = Chem.MolFromSmiles(product)
        fp = np.array(
            AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits=fpsize, useFeatures=useFeatures))
        if pfp is None:
            pfp = fp
        else:
            pfp += fp

    if pfp is not None and rfp is not None:
        pfp -= rfp
    return pfp


def reaction_condition(db, r_id, max_int, max_context, singleSlvt=True, with_smiles=True):
    """Return a set of conditions from instances associated with the reaction

    :param db: database info
    :param r_id: _id for a reaction
    :param max_int: the maximum number of instances to check for a reaction
    :param max_context: the maximum number of unqiue contexts per reaction
    :param singleSlvt: whether only use the first solvent
    :return: a list of contexts as tuples (T, solvent, reagent, catalyst, rxn_time, rxn_yield)
    """
    REACTION_DB = db['reactions']
    INSTANCE_DB = db['instances']
    CHEMICAL_DB = db['chemicals']
    SOLVENT_DB = db['solvents']

    rxn_doc = REACTION_DB.find_one({'_id': int(r_id)})

    # Look for conditions in the corresponding instances
    rxd_id_list = ['{}-{}'.format(rxn_doc['_id'], j) for j in range(1, int(rxn_doc['RX_NVAR']) + 1)]
    context_set = set()
    contexts = []

    for i, rxd_id in enumerate(rxd_id_list):
        num_context = len(context_set)
        if i >= max_int and len(contexts) > 0:
            break

        inst_doc = INSTANCE_DB.find_one({'_id': rxd_id})
        if not inst_doc:
            continue

        context_info = ''
        # Temp
        T = string_or_range_to_float(inst_doc['RXD_T'])
        if T is None or T == -1:  # skip if T was unparseable or not recorded
            continue

        # Solvent(s)
        solvent = ''
        context_info += 'solv:'
        for xrn in inst_doc['RXD_SOLXRN']:
            slvt = CHEMICAL_DB.find_one({'_id': xrn})
            if not slvt:
                print('########## COULD NOT FIND SOLVENT {} ###########'.format(xrn))
                continue
            smi_or_name = str(slvt['SMILES'])
            if (not smi_or_name) and with_smiles:
                continue
            if not smi_or_name:
                smi_or_name = str(slvt['IDE_CN'])
            solvent += smi_or_name + '.'
            context_info += str(slvt['IDE_CN']) + '(' + str(slvt['SMILES']) + ')' + ','
            if singleSlvt:
                break
        if not solvent:
            continue

        # Reagents
        reagent = ''
        context_info += 'rgt:'
        for xrn in inst_doc['RXD_RGTXRN']:
            rgt = CHEMICAL_DB.find_one({'_id': xrn})
            if not rgt:
                print('########## COULD NOT FIND REAGENT {} ###########'.format(xrn))
                continue
            smi_or_name = str(rgt['SMILES'])
            if (not smi_or_name) and with_smiles:
                continue
            if not smi_or_name:
                smi_or_name = str(rgt['IDE_CN'])
            reagent += smi_or_name + '.'
            context_info += str(rgt['IDE_CN']) + '(' + str(rgt['SMILES']) + ')' + ','

        # Catalysts
        catalyst = ''
        context_info += 'cat:'
        for xrn in inst_doc['RXD_CATXRN']:
            cat = CHEMICAL_DB.find_one({'_id': xrn})
            if not cat:
                print('########## COULD NOT FIND CATALYST {} ###########'.format(xrn))
                continue
            smi_or_name = str(cat['SMILES'])
            if (not smi_or_name) and with_smiles:
                continue
            if not smi_or_name:
                smi_or_name = str(cat['IDE_CN'])
            catalyst += smi_or_name + '.'
            context_info += str(cat['IDE_CN']) + '(' + str(cat['SMILES']) + ')' + ','

        if solvent and solvent[-1] == '.':
            solvent = solvent[:-1]
        if reagent and reagent[-1] == '.':
            reagent = reagent[:-1]
        if catalyst and catalyst[-1] == '.':
            catalyst = catalyst[:-1]

        # Time, yield
        rxn_time = inst_doc['RXD_TIM']
        rxn_yield = inst_doc['RXD_NYD']

        context_set.add(context_info)
        if len(context_set) > num_context:
            contexts.append((T, solvent, reagent, catalyst, rxn_time, rxn_yield))
        # context_info += 'T:{}C'.format(T) + ',t:{}min'.format(rxn_time) + ',y:{}%'.format(rxn_yield)

        if len(context_set) >= max_context:
            break
    return contexts


def n_rxn_condition(n, db, max_total_context, dists, idx, rxn_ids, max_int, max_context,
                    singleSlvt=True, with_smiles=True, dist_limit=0.3):
    """Reaction condition list from the top 10 NN

    :param n: int, the number of nearest neighbors to extract rxn conditions from, n <= 10 here
    :param db: database info
    :param dists, idx: np.array output from NearestNeighbor model for one rxn (10L, )
    :param rxn_ids: the list of reaction IDs for all the reactions in the model
    :return: A list of reaction conditions [(T, solvent, reagent, catalyst, rxn_time, rxn_yield), ()]
    """
    if n > int(idx.shape[0]):
        print('More rxn condition options requested than the number of NN, n is set to {}'.format(idx.shape[1]))
    if dists[0] > dist_limit:
        print('No neighbor is found within a cosine distance of {}'.format(dist_limit))

    contexts = []
    # int_ids = []
    for i, rid in enumerate(idx):
        if i >= n:
            break
        context = reaction_condition(db, rxn_ids[rid], max_int, max_context, singleSlvt, with_smiles)
        [contexts.append(c) for c in context]
        if len(contexts) >= max_total_context:
            contexts = contexts[:max_total_context]
            break
        # int_ids.append(rxn_ids[rid])
    return contexts


class NNConditionPredictor():
    """Reaction condition predictor based on Nearest Neighbor method"""

    def __init__(self, singleSlvt=True, with_smiles=True):
        """
        :param singleSlvt:
        :param with_smiles:
        """
        self.nnModel = None
        self.rxn_ids = []
        self.num_cond = 1
        self.dist_limit = 0.3
        self.singleSlvt = singleSlvt
        self.with_smiles = with_smiles
        self.max_total_context = 10
        self.max_int = 15
        self.max_context = 2

    def load_db_model(self, db_client, model_dir):
        """Loads the db client and the nearest neighbor model"""

        # client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
        # db = client['reaxys']
        self.db = db_client

        # Load the nearest neighbor model
        with open(os.path.join(model_dir, 'fpNN-10_BT_256Full.pickle'), 'rb') as infile:
            self.nnModel = joblib.load(infile)

        # Load the rxn ids associated with the nearest neighbor model
        rxd_ids = []
        rxn_ids = []
        with open(os.path.join(model_dir, 'RxnID_infoFull.txt'), 'r') as infile:
            rxn_ids.append(infile.readlines()[1:])  # a list of str(rxn_ids) with '\n'
        for id in rxn_ids[0]:
            rxd_ids.append(id.replace('\n', ''))
        self.rxn_ids = rxd_ids

    def load_predictor(self, userInput):
        """Loads the predictor based on user input"""
        self.num_cond = userInput['num_cond']
        self.dist_limit = userInput['dist_limit']
        self.singleSlvt = userInput['first_solvent_only']
        self.with_smiles = userInput['with_smiles_only']
        self.max_total_context = userInput['max_total_context']
        self.max_int = userInput['max_int']
        self.max_context = userInput['max_context']

    def step_n_conditions(self, n, rxn):
        """Reaction condition recommendations for a rxn (SMILES) from top n NN"""
        try:
            rxn_fp = create_rxn_Morgan2FP(rxn, fpsize=256, useFeatures=True)
            dists, ids = self.nnModel.kneighbors(rxn_fp)  # (1L, 10L)
            # # print info of neighbors
            # print('No.\tRxn ID\tDist')
            # for i, dist in enumerate(dists[0]):
            #     print('{}\t{}\t{}'.format(i, self.rxn_ids[ids[0][i]], dist))
            return n_rxn_condition(n, db=self.db, max_total_context=self.max_total_context, dists=dists[0], idx=ids[0],
                                   rxn_ids=self.rxn_ids, max_int=self.max_int, max_context=self.max_context,
                                   singleSlvt=self.singleSlvt, with_smiles=self.with_smiles, dist_limit=self.dist_limit)
        except Exception as e:
            print('Rxn {} with an exception: {}'.format(rxn, e))
            return None

    def path_condition(self, n, path):
        """Reaction condition recommendation for a reaction path with multiple reactions

            path: a list of reaction SMILES for each step
            return: a list of reaction context with n options for each step
        """
        rxn_fps = []
        for rxn in path:
            rxn_fps.append(create_rxn_Morgan2FP(rxn, fpsize=256, useFeatures=True))
        rxn_fps = np.array(rxn_fps)
        dists, ids = self.nnModel.kneighbors(rxn_fps)
        contexts = []
        for i, dist in enumerate(dists):
            context = []
            try:
                context = n_rxn_condition(n, db=self.db, max_total_context=self.max_total_context, dists=dist,
                                          idx=ids[i], rxn_ids=self.rxn_ids, max_int=self.max_int,
                                          max_context=self.max_context, singleSlvt=self.singleSlvt,
                                          with_smiles=self.with_smiles, dist_limit=self.dist_limit)
                contexts.append(context)
            except Exception as e:
                print('Step {} with an exception: {}'.format(i, e))
            
        return contexts
