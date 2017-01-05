import networkx as nx
from networkx.algorithms import bipartite
from pymongo import MongoClient
from Queue import Queue
import rdkit.Chem as Chem
from collections import defaultdict

class NetworkSearcher():
	'''This class is meant to be used to search the scraped
	set of reactions from Reaxys to look for known literature pathways'''

	def __init__(self, depth = 1):

		self.G = nx.DiGraph()
		self.expanded_xrns = defaultdict(lambda: depth)
		self.depth = depth
		self.queue = Queue()
		self.target = ''

		client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
		db = client['reaxys']
		self.CHEMICAL_DB = db['chemicals']
		self.REACTION_DB = db['reactions']

	def build(self, smiles):
		'''For a given SMILES (found in Reaxys), this method builds up
		the network'''

		self.target = smiles

		# Add the target compound
		chem_doc = self.CHEMICAL_DB.find_one({'SMILES': smiles}, ['_id'])

		if not chem_doc: # try to canonicalize
			mol = Chem.MolFromSmiles(smiles)
			if not mol: return False
			smiles = Chem.MolToSmiles(mol)
			chem_doc = self.CHEMICAL_DB.find_one({'SMILES': smiles}, ['_id'])

			if not chem_doc: return False
		
		# Save canonical SMILES that worked
		self.target = smiles

		xrn = chem_doc['_id']
		print('Matched target to ID {}'.format(xrn))
		self.G.add_node(xrn, smiles = smiles, distance = 0, buyable = self.is_buyable(xrn), bipartite = 0)
		self.queue.put(xrn)

		# Build
		while not self.queue.empty():
			xrn = self.queue.get() # XRN
			print(self.G.node[xrn])
			this_distance = self.G.node[xrn]['distance']
			if this_distance == self.depth: continue

			# Look for all reactions where this is the product
			for rx_doc in self.REACTION_DB.find({'RX_PXRN': xrn}, ['RX_RXRN']):

				# Add reaction node
				new_xrns = rx_doc['RX_RXRN']
				rx_id_string = 'rx_{}'.format(rx_doc['_id'])
				self.G.add_node(rx_id_string, bipartite = 1)

				# Add chemical nodes
				for new_xrn in new_xrns:
					# Don't add node if it has already been added with equal or closer distance
					if new_xrn in self.G.node:
						if self.G.node[new_xrn]['distance'] <= this_distance + 1: 
							continue

					# Look up SMILES
					chem_doc = self.CHEMICAL_DB.find_one({'_id': new_xrn}, ['SMILES'])
					if not chem_doc: 
						print('Warning: could not find entry for ID {}'.format(new_xrn))
						continue
					self.G.add_node(new_xrn, smiles = chem_doc['SMILES'], distance = this_distance + 1, bipartite = 0, 
						buyable = self.is_buyable(xrn))

					if this_distance + 1 < self.expanded_xrns[new_xrn] and not self.G.node[new_xrn]['buyable']:
						self.queue.put(new_xrn)
						self.expanded_xrns[new_xrn] = this_distance + 1
						print('Added {} to queue with distance {}'.format(new_xrn, this_distance + 1))

				# Add edges
				for new_xrn in new_xrns:
					self.G.add_edge(new_xrn, rx_id_string)
				self.G.add_edge(rx_id_string, xrn)

		return True # no errors

	def is_buyable(self, xrn):
		'''Checks to see if a compound is buyable by XRN'''
		chem_doc = self.CHEMICAL_DB.find_one({'_id': xrn}, ['SMILES', 'buyable_id'])
		if 'buyable_id' in chem_doc:
			return True
		return False

	def draw(self, label = None):
		import matplotlib.pyplot as plt

		if label == None:
			label = self.target

		nx.draw(self.G)
		plt.savefig('graph_for_{}.png'.format(label))
		nx.write_gexf(self.G, 'graph_for_{}.gexf'.format(label))

if __name__ == '__main__':
	NS = NetworkSearcher(depth = 3)
	NS.build('N/C(N)=C(\\[N+](=O)[O-])[N+](=O)[O-]')
	NS.draw('FOX-7')