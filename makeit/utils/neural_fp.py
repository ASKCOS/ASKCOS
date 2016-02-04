# Utilities for neural fingerprints
# based on description in http://arxiv.org/pdf/1509.09292v2.pdfv

import numpy as np
import rdkit.Chem.AllChem as AllChem
att_dtype = np.float32

class Graph():
	'''Describes an undirected graph class'''
	def __init__(self):
		self.nodes = []
		self.edges = []
		return

class Node():
	'''Describes an attributed node in an undirected graph'''
	def __init__(self, i = None, attributes = np.array([], dtype = att_dtype)):
		self.i = i
		self.attributes = attributes
		self.neighbors = [] # (atom index, bond index)
		return

class Edge():
	'''Describes an attributed edge in an undirected graph'''
	def __init__(self, connects = (), i = None, attributes = np.array([], dtype = att_dtype)):
		self.i = i
		self.attributes = attributes
		self.connects = connects # (atom index, atom index)
		return

def molToGraph(rdmol):
	'''Converts an RDKit molecule to an attributed undirected graph'''
	# Initialize
	graph = Graph()
	# Add atoms
	for atom in rdmol.GetAtoms():
		node = Node()
		node.i = atom.GetIdx()
		node.attributes = atomAttributes(atom)
		for neighbor in atom.GetNeighbors():
			node.neighbors.append((
				neighbor.GetIdx(),
				rdmol.GetBondBetweenAtoms(
					atom.GetIdx(),
					neighbor.GetIdx()
				).GetIdx()
			))
		graph.nodes.append(node)
	# Add bonds
	for bond in rdmol.GetBonds():
		edge = Edge()
		edge.i = bond.GetIdx()
		edge.attributes = bondAttributes(bond)
		edge.connects = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
		graph.edges.append(edge)
	return graph 

def bondAttributes(bond):
	'''Returns a numpy array of attributes for an RDKit bond

	From Neural FP defaults:
	The bond features were a concatenation of whether the bond type was single, double, triple,
	or aromatic, whether the bond was conjugated, and whether the bond was part of a ring.
	'''
	# Initialize
	attributes = []
	# Add bond type
	attributes += oneHotVector(
		bond.GetBondTypeAsDouble(),
		[1.0, 1.5, 2.0, 3.0]
	)
	# Add if is aromatic
	attributes.append(bond.GetIsAromatic())
	# Add if bond is conjugated
	attributes.append(bond.GetIsConjugated())
	# Add if bond is part of ring
	attributes.append(bond.IsInRing())

	return np.array(attributes, dtype = att_dtype)

def atomAttributes(atom):
	'''Returns a numpy array of attributes for an RDKit atom

	From ECFP defaults:
	<IdentifierConfiguration>
        <Property Name="AtomicNumber" Value="1"/>
        <Property Name="HeavyNeighborCount" Value="1"/>
        <Property Name="HCount" Value="1"/>
        <Property Name="FormalCharge" Value="1"/>
        <Property Name="IsRingAtom" Value="1"/>
    </IdentifierConfiguration>
    '''
	# Initialize
	attributes = []
	# Add atomic number (todo: finish)
	attributes += oneHotVector(
		atom.GetAtomicNum(),
		[1, 4, 5, 6, 7, 8, 9, 10, 'other']
	)
	# Add heavy neighbor count
	attributes += oneHotVector(
		len(atom.GetNeighbors()),
		[0, 1, 2, 3, 4, 5, 6]
	)
	# Add hydrogen count
	attributes += oneHotVector(
		atom.GetTotalNumHs(),
		[0, 1, 2, 3, 4]
	)
	# Add formal charge
	attributes.append(atom.GetFormalCharge())
	# Add boolean if in ring
	attributes.append(atom.IsInRing())
	# Add boolean if aromatic atom
	attributes.append(atom.GetIsAromatic())

	return np.array(attributes, dtype = att_dtype)

def oneHotVector(val, lst):
	'''Converts a value to a one-hot vector based on options in lst'''
	if val not in lst:
		return lst[-1]
	return map(lambda x: x == val, lst)

m = AllChem.MolFromSmiles('CC(O)C=C')
g = molToGraph(m)