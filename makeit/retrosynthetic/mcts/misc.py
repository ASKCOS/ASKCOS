import cPickle as pickle 
from scipy.sparse import csr_matrix, vstack as sparse_vstack 
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from multiprocessing import Pool
from functools import partial
import numpy as np 
import gzip 
import os

def mol_to_fp(mol, FINGERPRINT_SIZE = 16384, FP_rad = 3):
    if mol is None:
        return np.zeros((FINGERPRINT_SIZE,), dtype=np.int)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, FP_rad, nBits=FINGERPRINT_SIZE,
                                        useChirality=True, useFeatures=False), dtype=np.int)

def get_feature_vec(smiles, FINGERPRINT_SIZE, FP_rad):
    if not smiles:
        return np.zeros((FINGERPRINT_SIZE,), dtype=np.int)
    return mol_to_fp(Chem.MolFromSmiles(smiles), FINGERPRINT_SIZE, FP_rad)

def arr_for_cost(cost, SCALE = 1.0):
	val = np.array((float(cost),))
	return val

def get_states(pair, FINGERPRINT_SIZE, FP_rad, SCALE = 1.0):
	_id, smiles, cost, is_buyable = pair
	_id = int(_id)
	fps = get_feature_vec(smiles, FINGERPRINT_SIZE, FP_rad)
	carr = arr_for_cost(cost, SCALE)
	return (smiles, _id, fps, carr, is_buyable)

def save_sparse_tree(smile, smile_id, array, value, success, fid, FP_size):   
	array = csr_matrix(np.array(map(int, array)).astype(bool), (1,FP_size), dtype=bool)     
	matrix_parameters = [array.data, array.indices, array.indptr, array.shape, smile, smile_id, value, success]
	pickle.dump(matrix_parameters,fid,pickle.HIGHEST_PROTOCOL)

def value_network_training_states(smiles_id, Chemicals, Reactions, FP_rad = 3, FPS_size = 16384, fileName = ""):
	# Chemical
	smiles, states, values = [], [], []
	for chem_key, chem_dict in Chemicals.items():
		cost = float(chem_dict.cost)
		if cost > 0.0 and cost < 1000.0:
			smi, depth = chem_key
			smiles.append(smi)
			states.append(get_feature_vec(smi,FPS_size,FP_rad))
			values.append(cost)
	with open(fileName, "a+b") as fid:
		for smile, state, value in zip(smiles,states,values):
			save_sparse_tree(smile, smiles_id, state, float(value), 1, fid, FPS_size)
	print "... saved {} buyable chemicals to file.".format(len(values))
	# Reaction

def network_statistics(smiles_id, Chemicals, Reactions):
	# Log network statistics ...
	stats_location = "crn/mol_{}.stats".format(smiles_id)
	with gzip.open(stats_location, "a+") as fid:
		for key, Chem in Chemicals.items():
			smi, dep = key 
			_cost = Chem.cost
			_paths = Chem.counter 
			_visits = Chem.visit_count
			_incoming = len(Chem.incoming_reactions)
			_outgoing = len(Chem.outgoing_reactions)
			pickle.dump([smi,smiles_id,dep,_cost,_paths,_visits,_incoming,_outgoing], fid, pickle.HIGHEST_PROTOCOL)
	with gzip.open(stats_location, "a+") as fid:
		for key, Reac in Reactions.items():
			smi, dep = key 
			_cost = Reac.cost
			_paths = Reac.counter
			_visits = Reac.visit_count
			_incoming = len(Reac.incoming_chemicals)
			_outgoing = len(Reac.outgoing_chemicals)
			pickle.dump([smi,smiles_id,dep,_cost,_paths,_visits,_incoming,_outgoing], fid, pickle.HIGHEST_PROTOCOL)

