import numpy as np
# import rdkit.Chem.Descriptors as Descriptors
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
import rdkit.Chem.EState as EState
import rdkit.Chem.rdPartialCharges as rdPartialCharges
import rdkit.Chem.rdChemReactions as rdRxns
from makeit.utilities.atoms import atom_dftb
from makeit.utilities.io.logging import MyLogger
import rdkit.Chem as Chem 
descriptors_loc = 'descriptors'

att_dtype = np.float32

def oneHotVector(val, lst):
	'''Converts a value to a one-hot vector based on options in lst'''
	if val not in lst:
		val = lst[-1]
	return list(map(lambda x: x == val, lst))

def rxn_level_descriptors(rxn):
	'''
	Given an RDKit reaction, returns a reaction fingerprint as a numpy array
	**for very basic testing**
	'''
	settings = rdRxns.ReactionFingerprintParams()
	return np.array(rdRxns.CreateStructuralFingerprintForReaction(rxn, settings))

# def mol_level_descriptors(mol):
# 	'''
# 	Given an RDKit mol, returns a list of molecule-level descriptors 
# 	and their names

# 	returns: (labels, attributes)
# 	'''
	
# 	labels = [label for (label, f) in Descriptors._descList]
# 	attributes = [f(mol) for (label, f) in Descriptors._descList]
	
# 	return (labels, attributes)

def atom_level_descriptors(mol, include = ['functional'], asOneHot = False, ORIGINAL_VERSION = False):
	'''
	Given an RDKit mol, returns an N_atom-long list of lists,
	each of which contains atom-level descriptors and their names

	returns: (label, attributes)
	'''

	attributes = [[] for i in mol.GetAtoms()]
	labels = []
	if 'functional' in include:

		[attributes[i].append(x[0]) \
			for (i, x) in enumerate(rdMolDescriptors._CalcCrippenContribs(mol))]
		labels.append('Crippen contribution to logp')

		[attributes[i].append(x[1]) \
			for (i, x) in enumerate(rdMolDescriptors._CalcCrippenContribs(mol))]
		labels.append('Crippen contribution to mr')

		[attributes[i].append(x) \
			for (i, x) in enumerate(rdMolDescriptors._CalcTPSAContribs(mol))]
		labels.append('TPSA contribution')

		[attributes[i].append(x) \
			for (i, x) in enumerate(rdMolDescriptors._CalcLabuteASAContribs(mol)[0])]
		labels.append('Labute ASA contribution')

		[attributes[i].append(x) \
			for (i, x) in enumerate(EState.EStateIndices(mol))]
		labels.append('EState Index')

		rdPartialCharges.ComputeGasteigerCharges(mol)
		[attributes[i].append(float(a.GetProp('_GasteigerCharge'))) \
			for (i, a) in enumerate(mol.GetAtoms())]
		labels.append('Gasteiger partial charge')

		# Gasteiger partial charges sometimes gives NaN
		for i in range(len(attributes)):
			if np.isnan(attributes[i][-1]):
				attributes[i][-1] = 0.0

		[attributes[i].append(float(a.GetProp('_GasteigerHCharge'))) \
			for (i, a) in enumerate(mol.GetAtoms())]
		labels.append('Gasteiger hydrogen partial charge')

		# Gasteiger partial charges sometimes gives NaN
		for i in range(len(attributes)):
			if np.isnan(attributes[i][-1]):
				attributes[i][-1] = 0.0
	
	if 'structural' in include:
		[attributes[i].extend(atom_structural(mol.GetAtomWithIdx(i), asOneHot = asOneHot, ORIGINAL_VERSION = ORIGINAL_VERSION)) \
			for i in range(len(attributes))]
		labels.append('--many structural--')

	if 'dftb' in include:
		try:
			dftb_atom_atts = atom_dftb(mol)
		except ValueError as e:# often, an invalid element
			print(e)
			dftb_atom_atts = [[0 for i in range(18)] for j in range(mol.GetNumAtoms())]
		except KeyError as e: 
			print(e)
			dftb_atom_atts = [[0 for i in range(18)] for j in range(mol.GetNumAtoms())]
		[attributes[i].extend(dftb_atom_atts[i]) for i in range(mol.GetNumAtoms())]
		labels.append('--many DFTB--')

	return (labels, attributes)


def bond_structural(bond, asOneHot = False, extraOne = False):
	'''
	Returns a numpy array of attributes for an RDKit bond
	- Bond type as double
	- If bond is aromatic
	- If bond is conjugated
	- If bond is in ring
	'''

	# Redefine oneHotVector function
	if not asOneHot: 
		oneHotVectorFunc = lambda x: x[0]
	else:
		oneHotVectorFunc = oneHotVector

	# Initialize
	attributes = []
	# Add bond type
	attributes += oneHotVectorFunc(
		bond.GetBondTypeAsDouble(),
		[1.0, 1.5, 2.0, 3.0]
	)
	# Add if is aromatic
	attributes.append(bond.GetIsAromatic())
	# Add if bond is conjugated
	attributes.append(bond.GetIsConjugated())
	# Add if bond is part of ring
	attributes.append(bond.IsInRing())

	# NEED THIS FOR TENSOR REPRESENTATION - 1 IF THERE IS A BOND
	if extraOne: attributes.append(1)

	return np.array(attributes, dtype = att_dtype)

def atom_structural(atom, asOneHot = False, ORIGINAL_VERSION = False):
	'''
	Returns a numpy array of attributes for an RDKit atom
	- atomic number
	- number of heavy neighbors
	- total number of hydrogen neighbors
	- formal charge
	- is in a ring
	- is aromatic
	'''
	# Redefine oneHotVector function
	if not asOneHot: 
		oneHotVectorFunc = lambda x: x[0]
	else:
		oneHotVectorFunc = oneHotVector


	# Initialize
	attributes = []

	if ORIGINAL_VERSION: # F_atom = 32 mode
		attributes += oneHotVectorFunc(
			atom.GetAtomicNum(), 
			[5, 6, 7, 8, 9, 15, 16, 17, 35, 53, 999]
		)
		attributes += oneHotVectorFunc(
			len(atom.GetNeighbors()),
			[0, 1, 2, 3, 4, 5]
		)
		attributes += oneHotVectorFunc(
			atom.GetTotalNumHs(),
			[0, 1, 2, 3, 4]
		)
		attributes.append(atom.GetFormalCharge())
		attributes.append(atom.IsInRing())
		attributes.append(atom.GetIsAromatic())
		return np.array(attributes, dtype = att_dtype)

	
	# Add atomic number (todo: finish)
	attributes += oneHotVectorFunc(
		atom.GetAtomicNum(), 
		[3, 5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 35, 53, 999]
	)
	# Add heavy neighbor count
	attributes += oneHotVectorFunc(
		len(atom.GetNeighbors()),
		[0, 1, 2, 3, 4, 5]
	)
	# Add hydrogen count
	attributes += oneHotVectorFunc(
		atom.GetTotalNumHs(),
		[0, 1, 2, 3, 4]
	)
	# Add formal charge
	attributes.append(atom.GetFormalCharge())
	# Add boolean if in ring
	attributes.append(atom.IsInRing())
	# Add boolean if aromatic atom
	attributes.append(atom.GetIsAromatic())
	# Adjacent to aromatic ring but not aromatic itself
	attributes.append(atom.GetIsAromatic() == False and any([neighbor.GetIsAromatic() for neighbor in atom.GetNeighbors()]))

	# Halogen
	attributes.append(atom.GetAtomicNum() in [9, 17, 35, 53, 85, 117])
	# Chalcogen
	attributes.append(atom.GetAtomicNum() in [8, 16, 34, 52, 84, 116])
	# Pnictogens
	attributes.append(atom.GetAtomicNum() in [7, 15, 33, 51, 83])
	# Alkali
	attributes.append(atom.GetAtomicNum() in [3, 11, 19, 37, 55, 87])
	# Alkaline earth
	attributes.append(atom.GetAtomicNum() in [4, 12, 20, 38, 56, 88])
	# Common metals
	attributes.append(atom.GetAtomicNum() in [13, 22, 24, 25, 26, 27, 28, 29, 30, 33, 42, 44, 45, 46, 47, 48, 49, 50, 78, 80, 82])


	return np.array(attributes, dtype = att_dtype)

def edit_vector_lengths():
	mol = Chem.MolFromSmiles('[C:1][C:2]')
	(a, _, b, _) = edits_to_vectors((['1'],[],[('1','2',1.0)],[]), mol)
	F_atom = len(a[0])
	F_bond = len(b[0])
	return {'atoms': F_atom, 'bonds': F_bond}

def edits_to_vectors(edits, mol, atom_desc_dict = {}, return_atom_desc_dict = False, ORIGINAL_VERSION = False, include = ['functional', 'structural']):
	'''
	Given a set of edits (h_lost, h_gain, bond_lost, and bond_gain) from summarize_reaction_outcome,
	this function returns a set of vectors describing those edits.
	'''

	if not atom_desc_dict:
		atom_descriptors = atom_level_descriptors(mol, include = include, asOneHot = True, ORIGINAL_VERSION = ORIGINAL_VERSION)[1]
		atom_desc_dict = {a.GetProp('molAtomMapNumber'): atom_descriptors[i] for (i, a) in enumerate(mol.GetAtoms()) if a.HasProp('molAtomMapNumber')}
	if return_atom_desc_dict:
		return atom_desc_dict

	# h_lost, h_gain, bond_lost, bond_gain = edits
	return (
		[atom_desc_dict[molAtomMapNumber] for molAtomMapNumber in edits[0]], 
		[atom_desc_dict[molAtomMapNumber] for molAtomMapNumber in edits[1]], 
		[
		atom_desc_dict[molAtomMapNumber1] + 
		oneHotVector(bondOrder, [1.0, 1.5, 2.0, 3.0]) + 
		atom_desc_dict[molAtomMapNumber2] \
			for (molAtomMapNumber1, molAtomMapNumber2, bondOrder) in edits[2]
		],
		[
		atom_desc_dict[molAtomMapNumber1] + 
		oneHotVector(bondOrder, [1.0, 1.5, 2.0, 3.0]) + 
		atom_desc_dict[molAtomMapNumber2] \
			for (molAtomMapNumber1, molAtomMapNumber2, bondOrder) in edits[3]
		]
	)

def edits_to_tensor(candidate_edits, reactants, atom_desc_dict):
	
    Nc = len(candidate_edits)
    F_atom = edit_vector_lengths()['atoms']
    F_bond = edit_vector_lengths()['bonds']
    
    if Nc == 0:
        MyLogger.print_and_log('No candidate products found at all...?', descriptors_loc, level= 2)
        return None
    Ne1 = 1
    Ne2 = 1
    Ne3 = 1
    Ne4 = 1
    for (candidate_smiles, edits) in candidate_edits:
        Ne1 = max(Ne1, len(edits[0]))
        Ne2 = max(Ne2, len(edits[1]))
        Ne3 = max(Ne3, len(edits[2]))
        Ne4 = max(Ne4, len(edits[3]))
    x_h_lost = np.zeros((1, Nc, Ne1, F_atom))
    x_h_gain = np.zeros((1, Nc, Ne2, F_atom))
    x_bond_lost = np.zeros((1, Nc, Ne3, F_bond))
    x_bond_gain = np.zeros((1, Nc, Ne4, F_bond))

    smiles = []
    for c, (candidate_smiles, edits) in enumerate(candidate_edits):
        smiles.append(candidate_smiles)
        edit_h_lost_vec, edit_h_gain_vec, \
            edit_bond_lost_vec, edit_bond_gain_vec = edits_to_vectors(edits, reactants, 
                atom_desc_dict=atom_desc_dict)

        for (e, edit_h_lost) in enumerate(edit_h_lost_vec):
            x_h_lost[0, c, e, :] = edit_h_lost
        for (e, edit_h_gain) in enumerate(edit_h_gain_vec):
            x_h_gain[0, c, e, :] = edit_h_gain
        for (e, edit_bond_lost) in enumerate(edit_bond_lost_vec):
            x_bond_lost[0, c, e, :] = edit_bond_lost
        for (e, edit_bond_gain) in enumerate(edit_bond_gain_vec):
            x_bond_gain[0, c, e, :] = edit_bond_gain

    # Get rid of NaNs/Infs
    x_h_lost[np.isnan(x_h_lost)] = 0.0
    x_h_gain[np.isnan(x_h_gain)] = 0.0
    x_bond_lost[np.isnan(x_bond_lost)] = 0.0
    x_bond_gain[np.isnan(x_bond_gain)] = 0.0
    x_h_lost[np.isinf(x_h_lost)] = 0.0
    x_h_gain[np.isinf(x_h_gain)] = 0.0
    x_bond_lost[np.isinf(x_bond_lost)] = 0.0
    x_bond_gain[np.isinf(x_bond_gain)] = 0.0
    return [x_h_lost, x_h_gain, x_bond_lost, x_bond_gain]
'''
def edits_to_tensor(edits, reactants, lenAtom = 43, lenBond = 90, atom_desc_dict = None):
	
	(edit_h_lost_vec, edit_h_gain_vec, edit_bond_lost_vec, edit_bond_gain_vec) = edits_to_vectors(edits, reactants, atom_desc_dict = atom_desc_dict)
	x_h_lost = np.zeros((1, 1, len(edit_h_lost_vec), lenAtom))
	x_h_gain = np.zeros((1, 1, len(edit_h_gain_vec), lenAtom))
	x_bond_lost = np.zeros((1, 1, len(edit_bond_lost_vec), lenBond))
	x_bond_gain = np.zeros((1, 1, len(edit_bond_gain_vec), lenBond))

	for (e, edit_h_lost) in enumerate(edit_h_lost_vec):
		x_h_lost[0, 0, e, :] = edit_h_lost
	for (e, edit_h_gain) in enumerate(edit_h_gain_vec):
		x_h_gain[0, 0, e, :] = edit_h_gain
	for (e, edit_bond_lost) in enumerate(edit_bond_lost_vec):
		x_bond_lost[0, 0, e, :] = edit_bond_lost
	for (e, edit_bond_gain) in enumerate(edit_bond_gain_vec):
		x_bond_gain[0, 0, e, :] = edit_bond_gain

    # Get rid of NaNs
	x_h_lost[np.isnan(x_h_lost)] = 0.0
	x_h_gain[np.isnan(x_h_gain)] = 0.0
	x_bond_lost[np.isnan(x_bond_lost)] = 0.0
	x_bond_gain[np.isnan(x_bond_gain)] = 0.0
	x_h_lost[np.isinf(x_h_lost)] = 0.0
	x_h_gain[np.isinf(x_h_gain)] = 0.0
	x_bond_lost[np.isinf(x_bond_lost)] = 0.0
	x_bond_gain[np.isinf(x_bond_gain)] = 0.0
	
	return [x_h_lost, x_h_gain, x_bond_lost, x_bond_gain]
'''
