from makeit.utilities.fastfilter_utilities import Highway_self,pos_ct,true_pos, real_pos, set_keras_backend
from makeit.utilities.fingerprinting import create_rxn_Morgan2FP_separately
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from makeit.interfaces.scorer import Scorer
import h5py
import numpy as np
from scipy import sparse
import csv
import pandas as pd
from pymongo import MongoClient
from tqdm import tqdm
from keras.models import load_model
from keras import backend as K
import itertools, operator, random

scscore_prioritizer_loc = 'scscoreprioritizer'

# class Highway_self(Layer):

# 	def __init__(self, activation = 'elu',**kwargs):
# 		super(Highway_self, self).__init__(**kwargs)
# 		self.activation = activations.get(activation)
# 		self.transform_actv = activations.get('sigmoid')
		

# 	def build(self, input_shape):
# 		#weights of the dense layer
# 		self.kernel = self.add_weight(name = 'kernel',
# 								 shape = (input_shape[1],input_shape[1]),
# 								 initializer ='glorot_uniform',
# 								 trainable = True)
# 		self.bias = self.add_weight(name = 'bias',
# 								 shape = (input_shape[1],),
# 								 initializer ='zeros',
# 								 trainable = True)
# 		self.kernel_T = self.add_weight(name = 'kernel_T',
# 								 shape = (input_shape[1],input_shape[1]),
# 								 initializer ='glorot_uniform',
# 								 trainable = True)
# 		self.bias_T = self.add_weight(name = 'bias_T',
# 								 shape = (input_shape[1],),
# 								 initializer ='zeros',
# 								 trainable = True)
# 		self.input_dim = input_shape[1]
# 		# print(self.input_dim)
# 		super(Highway_self, self).build(input_shape)
	
# 	def call(self, x):
# 		transform_fun = self.activation(K.bias_add(K.dot(x,self.kernel), self.bias))
# 		transform_gate = self.transform_actv(K.bias_add(K.dot(x,self.kernel_T), self.bias_T))
# 		carry_gate = K.ones(self.input_dim,) - transform_gate
# 		output = transform_fun*transform_gate + x*carry_gate
# 		return output

# 	def compute_output_shape(self, input_shape):
# 		return (input_shape[0],input_shape[1])

# def pos_ct(y_true, y_pred):
# 	pos_pred = K.sum(gt((K.clip(y_pred, 0, 1)),0.5))
# 	return pos_pred
# def true_pos(y_true, y_pred):
# 	true_pos_ct = K.sum(gt((K.clip(y_pred*y_true, 0, 1)),0.5))
# 	return true_pos_ct
# def real_pos(y_true, y_pred):
# 	real_pos_ct = K.sum(gt((K.clip(y_true, 0, 1)),0.5))
# 	return real_pos_ct

# def create_rxn_Morgan2FP(rsmi, psmi, rxnfpsize = 2048, pfpsize=2048, useFeatures = False,calculate_rfp = True):
#     """Create a rxn Morgan (r=2) fingerprint as bit vector from SMILES string lists of reactants and products"""
#     # Modified from Schneider's code (2014)
#     if calculate_rfp is True:
# 	    rsmi = rsmi.encode('utf-8')
# 	    try:
# 	    	mol = Chem.MolFromSmiles(rsmi)
# 	    except Exception as e:
	    	
# 	    	return
# 	    # if mol is None:
# 	    # 	print(react)
# 	    try:
# 	        fp_bit = AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits = rxnfpsize, useFeatures=useFeatures)
# 	        fp = np.empty(rxnfpsize,dtype = 'int8')
# 	        DataStructs.ConvertToNumpyArray(fp_bit,fp)
# 	        # print(fp.dtype)
# 	        # fp = np.asarray(fp_bit)
# 	        # fp = AllChem.GetMorganFingerprint(mol=mol, radius=2, useFeatures=useFeatures)

# 	    except Exception as e:
# 	        print("Cannot build reactant fp due to {}".format(e))

# 	        return
	        
# 	    rfp = fp
#     else:
# 	    rfp = None

#     psmi = psmi.encode('utf-8')
#     try:
#     	mol = Chem.MolFromSmiles(psmi)
#     except Exception as e:
#     	print(psmi)
#     	return
#     # if mol is None:
#     # 	print(product)
#     try:
#         fp_bit = AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits = pfpsize, useFeatures=useFeatures)
#         fp = np.empty(pfpsize,dtype = 'int8')
#         DataStructs.ConvertToNumpyArray(fp_bit,fp)
#         # fp = np.asarray(fp_bit)
#         # fp = AllChem.GetMorganFingerprint(mol=mol, radius=2, useFeatures=useFeatures)

#     except Exception as e:
#     	print("Cannot build product fp due to {}".format(e))
#     	return
        
#     pfp = fp
#     # pfp_for_rxn = pfp
#     # for product in psmi:
#     # 	product = product.encode('utf-8')
#     #     mol = Chem.MolFromSmiles(product)
#     #     if mol is None:
#     #     	print(product)
#     #     try:
#     #         fp = np.array(
#     #             AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits=rxnfpsize, useFeatures=useFeatures))
#     #     except Exception as e:
#     #         print("Cannot build product fp due to {}".format(e))
#     #     if pfp_for_rxn is None:
#     #         pfp_for_rxn = fp
#     #     else:
#     #         pfp_for_rxn += fp
#     # if pfp_for_rxn is not None and rfp is not None:
#     #     rxnfp = pfp_for_rxn - rfp
#     return [pfp, rfp]


class FastFilterScorer(Scorer):
	def __init__(self):
		self.model = None

	def set_keras_backend(self,backend):
		if K.backend() != backend:
			os.environ['KERAS_BACKEND'] = backend
			reload(K)
			assert K.backend() == backend

	def load(self, model_path):
		##load model
		self.model = load_model(model_path, custom_objects = {'Highway_self':Highway_self,'pos_ct':pos_ct, 'true_pos':true_pos, 'real_pos':real_pos})

		# file=h5py.File(weights_path,'r')
		# weight = []
		# for i in range(len(file.keys())):
		#     weight.append(file['weight'+str(i)][:])
		# self.model.set_weights(weight)

	def evaluate(self, reactant_smiles, target, **kwargs):
		[pfp,rfp] = create_rxn_Morgan2FP_separately(reactant_smiles,target,rxnfpsize = 2048, pfpsize = 2048, useFeatures = False)
		pfp = np.asarray(pfp, dtype = 'float32')		
		rfp = np.asarray(rfp, dtype = 'float32')
		rxnfp = pfp - rfp

		score = self.model.predict([pfp.reshape(1,2048),rxnfp.reshape(1,2048)])
		outcome = {'smiles':target, 
				'template_ids':[], 
				'num_examples':0
				}
		all_outcomes = []
		all_outcomes.append([{'rank': 1.0,
							'outcome': outcome,
							'score': score[0][0],
							'prob': score[0][0],
							}])
		return all_outcomes
	def filter_with_threshold(self, reactant_smiles, target, threshold):
		[pfp,rfp] = create_rxn_Morgan2FP_separately(reactant_smiles,target,rxnfpsize = 2048, pfpsize = 2048, useFeatures = False)
		pfp = np.asarray(pfp, dtype = 'float32')		
		rfp = np.asarray(rfp, dtype = 'float32')
		rxnfp = pfp - rfp

		score = self.model.predict([pfp.reshape(1,2048),rxnfp.reshape(1,2048)])
		filter_flag = (score > threshold)
		return filter_flag



if __name__ == "__main__":

	# PFP3_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_3_69k_prediction/pfp_dataset_2.npz"
	# RFP3_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_3_69k_prediction/rfp_dataset_2.npz"
	# pfp_dataset = sparse.load_npz(PFP3_FPATH)
	# rfp_dataset = sparse.load_npz(RFP3_FPATH)
	ff = FastFilterScorer()
	ff.load(model_path ="/home/hanyug/Make-It/makeit/rxn_bin_class_NN/my_model.h5")
	score = ff.evaluate('CCO.CC(=O)O','CCOC(=O)C')
	print(score)
	score = ff.evaluate('[CH3:1][C:2](=[O:3])[O:4][CH:5]1[CH:6]([O:7][C:8]([CH3:9])=[O:10])[CH:11]([CH2:12][O:13][C:14]([CH3:15])=[O:16])[O:17][CH:18]([O:19][CH2:20][CH2:21][CH2:22][CH2:23][CH2:24][CH2:25][CH2:26][CH2:27][CH2:28][CH3:29])[CH:30]1[O:31][C:32]([CH3:33])=[O:34].[CH3:35][O-:36].[CH3:38][OH:39].[Na+:37]','CCCCCCCCCCOC1OC(CO)C(O)C(O)C1O')
	print(score)
	score = ff.evaluate('CNC.Cc1ccc(S(=O)(=O)OCCOC(c2ccccc2)c2ccccc2)cc1','CN(C)CCOC(c1ccccc1)c2ccccc2')
	print(score)
	# for i in range(10):
	# 	pfp = pfp_dataset[i,:].todense()
	# 	rfp = rfp_dataset[i,:].todense()
	# 	pfp = np.asarray(pfp, dtype = 'float32')		
	# 	rfp = np.asarray(rfp, dtype = 'float32')
	# 	rxnfp = pfp - rfp
	# 	score = ff.model.predict([pfp,rxnfp])
	# 	print(score)


	RID1_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_1_200k_prediction/rid_1.pickle"
	BIN1_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_1_200k_prediction/outcome_dataset_1.pickle"
	RID2_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_2_100k_prediction/rid_2.pickle"
	BIN2_FPATH = "/home/hanyug/Make-It/makeit/rxn_bin_class_NN/multiprocessing/part_2_100k_prediction/outcome_dataset_2.pickle"
	import pickle
	with open(RID1_FPATH,'r') as RID1:
		rid_list_1 = pickle.load(RID1)
	with open(RID2_FPATH,'r') as RID2:
		rid_list_2 = pickle.load(RID2)
	#sample the test set
	rid_list_test = rid_list_1[int(0.9*len(rid_list_1)):] +  rid_list_2[int(0.9*len(rid_list_2)):]
	rid_list_test = list(set(rid_list_test))
	rid_sample = np.random.choice(rid_list_test, size = 100)
	client = MongoClient('mongodb://guest:guest@askcos2.mit.edu/admin', 27017)
	db = client['prediction']
	reaxys_db = client['reaxys_v2']
	candidate_db = db['reaxys_edits_v2']
	instance_db = reaxys_db['instances']
	reaction_db = reaxys_db['reactions']
	chemical_db = reaxys_db['chemicals']
	MINIMUM_MAXPUB_YEAR = 1940
	results_df = pd.DataFrame(columns = ['rid','smiles','true_rxn','score'])
	ctr = 0
	for rxn_id in tqdm(rid_sample):
		# index += 1
		doc = candidate_db.find_one({'_id':str(rxn_id)+'-1'})
		# print(rxn_id)
		try:
			instance = instance_db.find_one({"_id":str(rxn_id)+'-1'},['RXD_STP','RXD_RGTXRN'])
		except:
			print('cannot find instance, skip')
			continue
		# print(instance)
		### keep only single step reactions
		if not instance: continue
		if instance['RXD_STP'] != ['1']: 
			print('multistep rxn, skip')
			continue
		reaction = reaction_db.find_one({"_id":rxn_id},['RX_MAXPUB'])
		if not reaction: 
			print("cannot find reaction, skip")
			continue
		if reaction['RX_MAXPUB'] == -1: 
			print("no publication year specified, skip")
			continue
		if int(reaction['RX_MAXPUB'][0]) < MINIMUM_MAXPUB_YEAR:
			print("outdated reaction, skip")
			continue
		
		# print(reaction)
		try:
			rsmi = doc["reactant_smiles"]
		
		# print(rsmi)
		# raw_input("press any key to continue")

			psmi = doc["product_smiles_true"]
		except:
			continue
		if '.' in psmi:
			print("multiproduct, skip")
			continue
		
		try:
			rct_mol = Chem.MolFromSmiles(rsmi)		
			[atom.ClearProp('molAtomMapNumber')for \
					atom in rct_mol.GetAtoms() if atom.HasProp('molAtomMapNumber')]
			rsmi = Chem.MolToSmiles(rct_mol)		
		except:
			continue
		# print(rsmi)
		# raw_input("press any key to continue")
		try:			
			rgt_smiles = [chemical_db.find_one({'_id':chem_id}, ["SMILES"])["SMILES"] for chem_id in instance['RXD_RGTXRN']]
		except:
			continue
		rgt_smiles = list(itertools.chain.from_iterable([rgt.split('.') for rgt in rgt_smiles]))
		# print(rgt_smiles)
		# raw_input("press any key to continue")
		rsmi_split = rsmi.split('.')
		rsmi_no_rgt = [reactant for reactant in rsmi_split if reactant not in rgt_smiles]
		rsmi = '.'.join(rsmi_no_rgt)
		### restric to single product reactions
		# print(rsmi)
		# raw_input("press any key to continue")

		csmis = doc["edit_candidates"]

		score = ff.evaluate(rsmi,psmi)
		results_df.loc[ctr] = [rxn_id, rsmi+'>>'+psmi, 1, score[0][0]["score"]]
		ctr+=1
		cdt_smiles = [candidate[0] for candidate in csmis \
				if candidate[0] != psmi]
		for csmi in cdt_smiles:
			score = ff.evaluate(rsmi,csmi)
			results_df.loc[ctr] = [rxn_id, rsmi+'>>'+csmi, 0, score[0][0]["score"]]
			ctr+=1

	results_df.to_csv('/home/hanyug/Make-It/makeit/rxn_bin_class_NN/results/example_results.csv')







	# with open('/home/hanyug/Make-It/makeit/NN_Distance/probs.csv') as data_src:
	# 	reader = csv.reader(data_src)
	# 	data_set  = list(reader)

	# ##delete title row
	# data_set = data_set[1:]

	# #create empty training and test set
	# rxn_smiles_list=[]
	# OPS_list=[]
	# RPS_list=[]

	# client = MongoClient('mongodb://guest:guest@askcos2.mit.edu/admin', 27017)

	# reaxys_db = client['reaxys_v2']
	# instance_db = reaxys_db['instances']
	# reaction_db = reaxys_db['reactions']
	# chemical_db = reaxys_db['chemicals']
	# i=0
	# for entry in tqdm(data_set[:1000]):
	# 	# print entry
	# 	# print entry[9]
	# 	# i+=1
	# 	# print(i)
	# 	##calculate Morgan2 finger prints for each 
	# 	# try:
	# 	rxn_id = int(entry[9].split('-')[0])
	# 	# print rxn_id
	# 	try:
	# 		rxn_smiles = reaction_db.find_one({'_id':rxn_id},['RXN_SMILES'])['RXN_SMILES']
	# 	except:
	# 		continue
	# 	splitted_rct_prd = rxn_smiles.split('>>')
	# 	rct = Chem.MolFromSmiles(splitted_rct_prd[0])
	# 	prd = Chem.MolFromSmiles(splitted_rct_prd[1])

	# 	# print(splitted_rct_prd)
	# 	[atom.ClearProp('molAtomMapNumber') for atom in rct.GetAtoms() if atom.HasProp('molAtomMapNumber')]
	# 	[atom.ClearProp('molAtomMapNumber') for atom in prd.GetAtoms() if atom.HasProp('molAtomMapNumber')]
	# 	rct = Chem.MolToSmiles(rct)
	# 	prd = Chem.MolToSmiles(prd)
	# 	# print(rct,prd)
	# 	if rct == None or prd == None:
	# 		continue

	# 	RPS = ff.evaluate(rct,prd)

	# 	rxn_smiles_list.append(entry[0])
	# 	OPS_list.append(entry[3])
	# 	RPS_list.append(RPS)

	# 	# except:
	# # 	# 	continue
	# # print(rxn_smiles_list)
	# # print(OPS_list)
	# # print(RPS_list)
	# RPS_mean = sum(RPS_list)/len(RPS_list)
	# print(RPS_mean)
	# df = pd.DataFrame({
	# 					'rxn_SMILES':rxn_smiles_list,
	# 					'Forward_predictor_score':OPS_list,
	# 					'Binary_class_score':RPS_list
	# 	})

	# df.to_csv('compare_bin_fwd_results.csv')