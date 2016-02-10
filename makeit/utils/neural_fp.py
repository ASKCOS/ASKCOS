# Utilities for neural fingerprints
# based on description in http://arxiv.org/pdf/1509.09292v2.pdfv

import numpy as np
import keras.backend as K
import rdkit.Chem.AllChem as AllChem
import copy
from theano.gof.type import Generic
att_dtype = np.float32

class Graph():
	'''Describes an undirected graph class'''
	def __init__(self):
		self.nodes = []
		self.num_nodes = 0
		self.edges = []
		self.num_edges = 0
		return

	def nodeAttributes(self):
		'''Returns 2D array where (#, :) contains attributes of node #'''
		return K.variable(np.vstack([x.attributes for x in self.nodes]))
	
	def edgeAttributes(self):
		'''Returns 2D array where (#, :) contains attributes of edge #'''
		return K.variable(np.vstack([x.attributes for x in self.edges]))

	def nodeNeighbors(self):
		return [x.neighbors for x in self.nodes]

	def clone(self):
		'''clone() method to trick Theano'''
		return copy.deepcopy(self) 

	def dump_as_tensor(self):
		'''Method to represent attributed graph as a giant tensor

		The tensor is N_node x N_node x N_attributes.

		For a given node, A_i,i,: is a vector of that node's features, followed
		  by zeros where edge attributes would be.
		For a pair of nodes i and j, A_i,j,: is a vector of node j's features, 
		  followed by the edge attributes connecting the two nodes.

		This representation is not as efficient as it could be, but we need to
		  pack the whole graph information into a single tensor in order to use
		  Keras/Theano easily'''

		# Bad input handling
		if not self.edges or not self.nodes:
			print('Error with graph!')
			print('edges: {}'.format(self.edges))
			print('nodes: {}'.format(self.nodes))

		N_nodes = len(self.nodes)
		N_features = sizeAttributeVector()
		tensor = np.zeros((N_nodes, N_nodes, N_features))
		edgeAttributes = np.vstack([x.attributes for x in self.edges])
		nodeAttributes = np.vstack([x.attributes for x in self.nodes])
		nodeNeighbors = self.nodeNeighbors()
		# Assign diagonal entries
		for i, node in enumerate(self.nodes):
			tensor[i, i, :] = np.concatenate((nodeAttributes[i], np.zeros_like(edgeAttributes[0])))
		# Assign bonds now
		for e, edge in enumerate(self.edges):
			(i, j) = edge.connects
			tensor[i, j, :] = np.concatenate((nodeAttributes[j], edgeAttributes[e]))
			tensor[j, i, :] = np.concatenate((nodeAttributes[i], edgeAttributes[e]))

		return tensor


class Node():
	'''Describes an attributed node in an undirected graph'''
	def __init__(self, i = None, attributes = np.array([], dtype = att_dtype)):
		self.i = i
		self.attributes = attributes # 1D array
		self.neighbors = [] # (atom index, bond index)
		return

class Edge():
	'''Describes an attributed edge in an undirected graph'''
	def __init__(self, connects = (), i = None, attributes = np.array([], dtype = att_dtype)):
		self.i = i
		self.attributes = attributes # 1D array
		self.connects = connects # (atom index, atom index)
		return

def molToGraph(rdmol):
	'''Converts an RDKit molecule to an attributed undirected graph'''
	# Initialize
	graph = Graph()
	# Add bonds
	for bond in rdmol.GetBonds():
		edge = Edge()
		edge.i = bond.GetIdx()
		edge.attributes = bondAttributes(bond)
		edge.connects = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
		graph.edges.append(edge)
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
	# Add counts, for convenience
	graph.num_edges = len(graph.edges)
	graph.num_nodes = len(graph.nodes)
	return graph 

def padGraphTensor(old_tensor, new_dsize):
	'''This function takes an input tensor of dsize x dsize x Nfeatures and pads 
	up the first two dimensions to new_dsize with zeros as needed'''
	
	old_shape = old_tensor.shape
	new_tensor = np.zeros((new_dsize, new_dsize, old_shape[2]))
	for i in range(old_shape[0]):
		for j in range(old_shape[1]):
			for k in range(old_shape[2]):
				new_tensor[i, j, k] = old_tensor[i, j, k]

	return new_tensor

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

	# NEED THIS FOR TENSOR REPRESENTATION - 1 IF THERE IS A BOND
	attributes.append(1)

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
		[5, 6, 7, 8, 9, 15, 16, 17, 35, 53, 999]
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
		val = lst[-1]
	return map(lambda x: x == val, lst)

def sizeAttributeVector():
	m = AllChem.MolFromSmiles('CC')
	g = molToGraph(m)
	a = g.nodes[0]
	b = g.edges[0]
	return len(a.attributes) + len(b.attributes)