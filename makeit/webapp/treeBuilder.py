from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import rdkit.Chem as Chem          
from rdkit.Chem import AllChem
import numpy as np
import Queue as VanillaQueue
from multiprocessing import Pool, cpu_count, Process, Manager, Queue
from functools import partial # used for passing args to multiprocessing
import time

def expansion_worker(i, expansion_queue, results_queue, paused, done):
	'''
	The target for a worker Process to expand to generate precursors

	Pulls (_id, smiles) from the expansion_queue
	Adds (_id, children) to the results_quere, 
	    where each child is a 2-tuple of (rxn_dict, mols_list)
	    and each mol in mols_list is just a smiles string
	'''
	while True:
		# If done, stop
		if done.value:
			print('Worker {} saw done signal, terminating'.format(i))
			break
		# If paused, wait and check again
		if paused.value:
			print('Worker {} saw pause signal, sleeping for 1 second'.format(i))
			time.sleep(1)
			continue
		# Grab something off the queue
		try:
			(_id, smiles) = expansion_queue.get()
			print('Worker grabbed {} (ID {}) to expand'.format(smiles, _id))

			# Fake processing
			children = []
			children.append((
				{'info': 'this is a reaction'},
				[
					smiles + '-a', 
					smiles + '-b'
				]
			))
			results_queue.put((_id, children))

			time.sleep(1)

		except VanillaQueue.Empty:
			pass


def coordinator(current_id, tree_dict, expansion_queue, results_queue, chem_to_id, paused, done):
	'''
	The coordinator
	'''
	while True:
		print('Coordinator done signal: {}'.format(done.value))
		# If done, stop
		if done.value:
			print('Coordinator saw done signal, stopping')
			break
		# If paused, wait and check again
		if paused.value:
			print('Coordinator saw pause signal, sleeping 1 second')
			time.sleep(1)
			continue
		try:
			(_id, children) = results_queue.get()
			print('Got results for node ID {}'.format(_id))

			# Assign unique number
			for (rxn, mols) in children:
				rxn_id = current_id.value
				current_id.value += 1 # this is only okay because there is ONE coordinator process

				# For the parent molecule, record child reactions
				tree_dict[_id]['is_product_of'].append(rxn_id)
			
				# For the reaction, keep track of children IDs
				tree_dict[rxn_id] = {
					'info': rxn['info'],
					'is_reaction_of': [],
				}

				for mol in mols:
					
					# New chemical?
					if mol not in chem_to_id:

						chem_id = current_id.value
						current_id.value += 1 # this is only okay because there is ONE coordinator process

						tree_dict[chem_id] = {
							'smiles': mol,
							'is_product_of': [],
						}
						chem_to_id[mol] = chem_id

						# Add to queue to get expanded
						expansion_queue.put((chem_id, mol))
						print('Put {} (ID {}) in expansion_queue'.format(mol, chem_id))

					else:
						chem_id = chem_to_id[mol]

					# Record
					tree_dict[rxn_id]['is_reaction_of'].append(chem_id)

		except VanillaQueue.Empty:
			pass


class TreeBuilder:
	'''
	The Transformer class defines an object which can be used to perform
	one-step retrosyntheses for a given molecule.
	'''

	def __init__(self, nb_workers = 2):

		# Number of workers to expand nodes
		self.nb_workers = nb_workers

	def start_building(self, smiles):
		'''
		Begin building the network
		'''

		# Queue of nodes that have yet to be expanded
		self.expansion_queue = Queue()
		# Queue of results that have yet to be processed and IDed
		self.results_queue = Queue()

		# Dictionary storing the overall tree
		self.manager = Manager()
		self.tree_dict = self.manager.dict()
		self.chem_to_id = self.manager.dict()
		self.paused = self.manager.Value('i', 0)
		self.done = self.manager.Value('i', 0)
		self.current_id = self.manager.Value('i', 1)
		self.workers = []

		# Initialize the queue
		self.current_id.value = 2
		self.expansion_queue.put((1, smiles))
		self.tree_dict[1] = {
			'smiles': smiles,
			'is_product_of': [],
		}

		# Begin processes
		for i in range(self.nb_workers):
			p = Process(target = expansion_worker, args = (i, self.expansion_queue, self.results_queue, self.paused, self.done))
			self.workers.append(p)
			p.start()

		self.coordinator = Process(target = coordinator, args = (self.current_id, self.tree_dict, self.expansion_queue, self.results_queue, self.chem_to_id, self.paused, self.done))
		self.coordinator.start()


	def stop_building(self):
		'''
		Stop building
		'''
		
		# Signal done
		self.done.value = 1

		print('Changed `done` signal to True')

		# Join up workers
		time.sleep(5)	
		for p in self.workers + [self.coordinator]:
			if p.is_alive():
				print('Process is still alive?')
				p.terminate()

		print('Stopped building, all processes done')

		return True