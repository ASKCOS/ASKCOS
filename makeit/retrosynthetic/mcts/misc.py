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
    """Converts a molecule to its fingerprint.

    Args:
        mol (Chem.rdchem.Mol): Molecule to get the fingerprint of.
        FINGERPRINT_SIZE (int, optional): Size of the fingerprint.
            (default: {16384})
        FP_rad (int, optional): Fingerprint radius. (default: {3})

    Returns:
        np.ndarray of np.int: Morgan fingerprint of the given molecule.
    """
    if mol is None:
        return np.zeros((FINGERPRINT_SIZE,), dtype=np.int)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, FP_rad, nBits=FINGERPRINT_SIZE,
                                        useChirality=True, useFeatures=False), dtype=np.int)

def get_feature_vec(smiles, FINGERPRINT_SIZE, FP_rad):
    """Converts a SMILES string to the fingerprint of the molecule.

    Args:
        smiles (str): SMILES string of molecule to get the fingerprint of.
        FINGERPRINT_SIZE (int): Size of the fingerprint.
        FP_rad (int): Fingerprint radius.

    Returns:
        np.ndarray of np.int: Morgan fingerprint of the given molecule.
    """
    if not smiles:
        return np.zeros((FINGERPRINT_SIZE,), dtype=np.int)
    return mol_to_fp(Chem.MolFromSmiles(smiles), FINGERPRINT_SIZE, FP_rad)

def arr_for_cost(cost, SCALE = 1.0):
    """Returns cost in an np.ndarray.

    Args:
        cost (float): Cost to be put into the array.
        SCALE (float, optional): Unused. (default: {1.0})
    """
    val = np.array((float(cost),))
    return val

def get_states(pair, FINGERPRINT_SIZE, FP_rad, SCALE = 1.0):
    """Gets information about a molecule.

    Args:
        pair (4-tuple of int, str, float, bool): ID, SMILES string, cost, and
            whether the molecule can be bought.
        FINGERPRINT_SIZE (int): Size of the fingerprint to return.
        FP_rad (int): Radius of the fingerprint.
        SCALE (float, optional): Passed to arr_for_cost; unused.
            (default: {1.0})

    Returns:
        5-tuple of (str, int, np.ndarray, np.ndarray, bool): SMILES string, ID,
            fingerprint, cost array, and whether the molecule is buyable.
    """
    _id, smiles, cost, is_buyable = pair
    _id = int(_id)
    fps = get_feature_vec(smiles, FINGERPRINT_SIZE, FP_rad)
    carr = arr_for_cost(cost, SCALE)
    return (smiles, _id, fps, carr, is_buyable)

def save_sparse_tree(smile, smile_id, array, value, success, fid, FP_size):
    """Saves a sparse tree for the molecule to disk.

    Args:
        smile (str): SMILES string of the given molecule.
        smile_id (??): ID for the molecule.
        array (np.ndarray): Fingerprint of molecule.
        value (float): Chemical cost.
        success (): ??
        fid (File Descriptor): File descriptor to save to.
        FP_size (int): Fingerprint size.
    """
    array = csr_matrix(np.array(map(int, array)).astype(bool), (1,FP_size), dtype=bool)
    matrix_parameters = [array.data, array.indices, array.indptr, array.shape, smile, smile_id, value, success]
    pickle.dump(matrix_parameters,fid,pickle.HIGHEST_PROTOCOL)

def value_network_training_states(smiles_id, Chemicals, Reactions, FP_rad = 3, FPS_size = 16384, fileName = ""):
    """Saves information about buyable chemicals to disk.

    Args:
        smiles_id (): ID for the molecule??
        Chemicals (dict of {(str, int): Chemical}): Chemcials to be saved.
        Reactions (dict): Chemcials to be saved. Unused.
        FP_rad (int, optional): Fingerprint radius. (default: {3})
        FPS_size (int, optional): Fingerprint size. (default: {16384})
        filename (str, optional): Filename to save to. (default: {''})
    """
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
    print("... saved {} buyable chemicals to file.".format(len(values)))
    # Reaction

def network_statistics(smiles_id, Chemicals, Reactions):
    """Logs network statistics.

    Args:
        smiles_id (): ??
        Chemicals (dict of {(str, int): Chemical}): Chemicals to store
            information about.
        Reactions (dict of {(str, int): Reaction}): Reactions to store
            information about.
    """
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
