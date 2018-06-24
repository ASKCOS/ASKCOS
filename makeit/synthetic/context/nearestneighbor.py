import makeit.global_config as gc
from pymongo import MongoClient
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import makeit.utilities.io.pickle as pickle
import numpy as np
import random
from collections import deque
# from sklearn.neighbors import NearestNeighbors as NN
from sklearn.externals import joblib
import makeit.utilities.strings as strings
import makeit.utilities.fingerprinting as fp
from makeit.utilities.io.logger import MyLogger
from makeit.interfaces.context_recommender import ContextRecommender
contextRecommender_loc = 'contextRecommender'

class NNContextRecommender(ContextRecommender):
    """Reaction condition predictor based on Nearest Neighbor method"""

    def __init__(self, max_contexts = 10, singleSlvt=True, with_smiles=True, done = None, REACTION_DB = None,
                 CHEMICAL_DB=None, INSTANCE_DB= None, cachelength=100):
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
        self.max_total_context = max_contexts
        self.max_int = 15
        self.max_context = 2
        self.done = done
        self.REACTION_DB = REACTION_DB
        self.CHEMICAL_DB = CHEMICAL_DB
        self.INSTANCE_DB = INSTANCE_DB
        self.cache = {}
        self.cache_q = deque([], maxlen=cachelength)
    
    def load(self, model_path = "", info_path = ""):
        self.load_databases()
        self.load_nn_model(model_path, info_path)
        
        MyLogger.print_and_log('Context recommender has been loaded.', contextRecommender_loc)
        
        #multiprocessing notify done
        if self.done == None:
            pass
        else:
            self.done.value = 1
    #Model path should be relative!
    def load_nn_model(self, model_path = "", info_path = ""):
        """Loads the nearest neighbor model"""
        
        
        if not model_path:      
            MyLogger.print_and_log('Cannot load nearest neighbor context recommender without a specific path to the model. Exiting...', contextRecommender_loc, level = 3)
        if not info_path:
            MyLogger.print_and_log('Cannot load nearest neighbor context recommender without a specific path to the model info. Exiting...', contextRecommender_loc, level = 3)
        # Load the nearest neighbor model
        
        with open(model_path, 'rb') as infile:
            self.nnModel = joblib.load(infile)
            
        # Load the rxn ids associated with the nearest neighbor model
        rxd_ids = []
        rxn_ids = []
        with open(info_path, 'r') as infile:
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

    def load_databases(self):
        db_client = MongoClient(gc.MONGO['path'], gc.MONGO['id'], connect = gc.MONGO['connect'])
        
        db = db_client[gc.REACTIONS['database']]
        self.REACTION_DB = db[gc.REACTIONS['collection']]
        
        db = db_client[gc.INSTANCES['database']]
        self.INSTANCE_DB = db[gc.INSTANCES['collection']]
        
        db = db_client[gc.CHEMICALS['database']]
        self.CHEMICAL_DB = db[gc.CHEMICALS['collection']]
        
    
    def get_n_conditions(self, rxn, n=10, singleSlvt=True, with_smiles=True):
        """
        Reaction condition recommendations for a rxn (SMILES) from top n NN
        Returns the top n parseable conditions.
        """
        self.singleSlvt = singleSlvt
        self.with_smiles = with_smiles

        #If the databases haven't been loaded into the object yet, do so.
        if not (self.REACTION_DB and self.INSTANCE_DB and self.CHEMICAL_DB):
            self.load_databases()
    
        try:
            if rxn in self.cache:
                return self.cache[rxn][:n]

            rxn_fp = fp.create_rxn_Morgan2FP(rxn, fpsize=256, useFeatures=True)
            dists, ids = self.nnModel.kneighbors(rxn_fp)  # (1L, 10L)
            # # print info of neighbors
            # print('No.\tRxn ID\tDist')
            # for i, dist in enumerate(dists[0]):
            #     print('{}\t{}\t{}'.format(i, self.rxn_ids[ids[0][i]], dist))
            
            # Grab as many contexts as allowed for caching
            conditions = self.n_rxn_condition(self.max_total_context, dists=dists[0], idx=ids[0])

            # Save in cache
            if len(self.cache_q) == self.cache_q.maxlen:
                rxn_to_remove = self.cache_q.pop()
                del self.cache[rxn_to_remove]
            self.cache_q.appendleft(rxn)
            self.cache[rxn] = conditions
            return conditions[:n]


        except Exception as e:

            MyLogger.print_and_log('Failed for reaction {} because {}. Returning None.'.format(rxn,e), contextRecommender_loc, level=2)
            
            return [[]]

    def path_condition(self, n, path):
        """Reaction condition recommendation for a reaction path with multiple reactions

            path: a list of reaction SMILES for each step
            return: a list of reaction context with n options for each step
        """
        rxn_fps = []
        for rxn in path:
            rxn_fps.append(fp.create_rxn_Morgan2FP(rxn, fpsize=256, useFeatures=True))
        rxn_fps = np.array(rxn_fps)
        dists, ids = self.nnModel.kneighbors(rxn_fps)
        contexts = []
        for i, dist in enumerate(dists):
            context = []
            try:
                context = self.n_rxn_condition(n, dists=dist, idx=ids[i])
                contexts.append(context)
            except Exception as e:
                print('Step {} with an exception: {}'.format(i, e))
            
        return contexts
    
    def n_rxn_condition(self, n, dists, idx):
        """
        Reaction condition list from the top 10 NN
        Only returns parseable solvents
        :param n: int, the number of nearest neighbors to extract rxn conditions from, n <= 10 here
        :param dists, idx: np.array output from NearestNeighbor model for one rxn (10L, )
        :return: A list of reaction conditions [(T, solvent, reagent, catalyst, rxn_time, rxn_yield), ()]
        """
        if n > int(idx.shape[0]):
            print('More rxn condition options requested than the number of NN, n is set to {}'.format(idx.shape[1]))
            if dists[0] > self.dist_limit:
                print('No neighbor is found within a cosine distance of {}'.format(self.dist_limit))
                
        contexts = []
        # int_ids = []
        for i, rid in enumerate(idx):
            if i >= n:
                break
            context = self.reaction_condition(self.rxn_ids[rid])
            for c in context:
                if not (c == []):
                    contexts.append(c)
            if len(contexts) >= self.max_total_context:
                contexts = contexts[:self.max_total_context]
                break
            # int_ids.append(self.rxn_ids[rid])
        
        return contexts
    
    def reaction_condition(self, r_id):
        """Return a set of conditions from instances associated with the reaction
    
        :param db: database info
        :param r_id: _id for a reaction
        :param max_int: the maximum number of instances to check for a reaction
        :param max_context: the maximum number of unqiue contexts per reaction
        :param singleSlvt: whether only use the first solvent
        :return: a list of contexts as tuples (T, solvent, reagent, catalyst, rxn_time, rxn_yield)
        """
        
        rxn_doc = self.REACTION_DB.find_one({'_id': int(r_id)})
        if not rxn_doc:
            return [[]]
        # Look for conditions in the corresponding instances
        rxd_id_list = ['{}-{}'.format(rxn_doc['_id'], j) for j in range(1, int(rxn_doc['RX_NVAR']) + 1)]
        context_set = set()
        contexts = []
        for i, rxd_id in enumerate(rxd_id_list):
            num_context = len(context_set)
            if i >= self.max_int and len(contexts) > 0:
                break
            
            inst_doc = self.INSTANCE_DB.find_one({'_id': rxd_id})
            if not inst_doc:
                continue
            
            context_info = ''
            # Temp
            T = strings.string_or_range_to_float(inst_doc['RXD_T'])
            if T is None or T == -1:  # skip if T was unparseable or not recorded
                continue
            
            # Solvent(s)
            solvent = ''
            context_info += 'solv:'
            for xrn in inst_doc['RXD_SOLXRN']:
                slvt = self.CHEMICAL_DB.find_one({'_id': xrn})
                if not slvt:
                    MyLogger.print_and_log('COULD NOT FIND SOLVENT {}'.format(xrn), contextRecommender_loc)
                    continue
                smior_name = str(slvt['SMILES'])
                if (not smior_name) and self.with_smiles:
                    continue
                if not smior_name:
                    smior_name = str(slvt['IDE_CN'])
                    solvent += smior_name + '.'
                    context_info += str(slvt['IDE_CN']) + '(' + str(slvt['SMILES']) + ')' + ','
                    if self.singleSlvt:
                        break
                    if not solvent:
                        continue

            # Reagents
            reagent = ''
            context_info += 'rgt:'
            for xrn in inst_doc['RXD_RGTXRN']:
                rgt = self.CHEMICAL_DB.find_one({'_id': xrn})
                if not rgt:
                    print('########## COULD NOT FIND REAGENT {} ###########'.format(xrn))
                    continue
                smior_name = str(rgt['SMILES'])
                if (not smior_name) and self.with_smiles:
                    continue
                if not smior_name:
                    smior_name = str(rgt['IDE_CN'])
                reagent += smior_name + '.'
                context_info += str(rgt['IDE_CN']) + '(' + str(rgt['SMILES']) + ')' + ','

            # Catalysts
            catalyst = ''
            context_info += 'cat:'
            for xrn in inst_doc['RXD_CATXRN']:
                cat = self.CHEMICAL_DB.find_one({'_id': xrn})
                if not cat:
                    print('########## COULD NOT FIND CATALYST {} ###########'.format(xrn))
                    continue
                smior_name = str(cat['SMILES'])
                if (not smior_name) and self.with_smiles:
                    continue
                if not smior_name:
                    smior_name = str(cat['IDE_CN'])
                catalyst += smior_name + '.'
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
           
            #Test that all proposed conditions are actually parseable!
            try:
                solvent_mols = [Chem.MolFromSmiles(sol) for sol in solvent.split('.')]
                if None in solvent_mols:
                    MyLogger.print_and_log('Unparseable solvent recommended. Not adding {}'.format(solvent),contextRecommender_loc)
                    continue
            except Exception as e:
                MyLogger.print_and_log('Unparseable solvent recommended. Not adding {}'.format(solvent),contextRecommender_loc)
                continue
            try:
                reagents_mols = [Chem.MolFromSmiles(rea) for rea in reagent.split('.')]
                if None in reagents_mols:
                    MyLogger.print_and_log('Unparseable reagent recommended. Not adding {}'.format(reagent), contextRecommender_loc)
                    continue
            except Exception as e:
                MyLogger.print_and_log('Unparseable reagent recommended. Not adding {}'.format(reagent), contextRecommender_loc)
                continue
            try:
                catalyst_mols = [Chem.MolFromSmiles(cat) for cat in catalyst.split('.')]
                if None in catalyst_mols:
                    MyLogger.print_and_log('Unparseable catalyst recommended. Not adding {}'.format(catalyst), contextRecommender_loc)
                    continue
            except Exception as e:
                MyLogger.print_and_log('Unparseable catalyst recommended. Not adding {}'.format(catalyst), contextRecommender_loc)
                continue
            context_set.add(context_info)
            if len(context_set) > num_context:
                contexts.append((T, solvent, reagent, catalyst, rxn_time, rxn_yield))
                # context_info += 'T:{}C'.format(T) + ',t:{}min'.format(rxn_time) + ',y:{}%'.format(rxn_yield)
                
            if len(context_set) >= self.max_context:
                break
        
        return contexts

if __name__ == '__main__':
    import global_config as gc
    cont = NNContextRecommender()
    cont.load_nn_model(model_path = gc.CONTEXT_REC['model_path'], info_path = gc.CONTEXT_REC['info_path'])
    print(cont.get_n_conditions('CN1C2CCC1CC(O)C2.O=C(O)C(CO)c1ccccc1>>CN1C2CCC1CC(C2)OC(=O)C(CO)c3ccccc3', 10))