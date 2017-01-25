from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import rdkit.Chem as Chem          
from rdkit.Chem import AllChem
import numpy as np
import Queue as VanillaQueue
from multiprocessing import Pool, cpu_count, Process, Manager, Queue
from functools import partial # used for passing args to multiprocessing
import time


def expansion_worker(i, expansion_queues, results_queue, paused, done, retrotransformer):
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
			#print('Worker {} saw pause signal, sleeping for 1 second'.format(i))
			time.sleep(1)
			continue
		# Grab something off the queue
		for j in range(len(expansion_queues))[::-1]:
			try:
				(_id, smiles) = expansion_queues[j].get(timeout = 0.05) # short timeout
				#print('Worker {} grabbed {} (ID {}) to expand from queue {}'.format(i, smiles, _id, j))


				results_queue.put((_id, get_children_from_smiles(retrotransformer, smiles)))
				#print('Worker {} added children of {} (ID {}) to results queue'.format(i, smiles, _id))
				
			except VanillaQueue.Empty:
				#print('Queue {} empty for worker {}'.format(j, i))
				pass


def get_children_from_smiles(retrotransformer, smiles):
	'''
	For a worker that just grabbed a smiles string, this function will apply templates
	and return a list of children 
	'''

	result = retrotransformer.perform_retro(smiles)
	precursors = result.return_top(n = 20)
	children = []
	for precursor in precursors:
		children.append((
			{
				'info': 'this is a reaction',
				'num_examples': precursor['num_examples'],
				'score': precursor['score'],
			},
			precursor['smiles_split']
		))

	return children


def coordinator(current_id, tree_dict, expansion_queues, results_queue, chem_to_id, paused, done, buyable_leaves, pricer):
	'''
	The coordinator

	TODO: add keeping track of buyable final nodes
	TODO: add multiple queueues
	'''
	while True:
		# If done, stop
		if done.value:
			print('Coordinator saw done signal, stopping')
			break
		# If paused, wait and check again
		if paused.value:
			#print('Coordinator saw pause signal, sleeping 1 second')
			time.sleep(1)
			continue
		try:
			(_id, children) = results_queue.get()
			#print('Got results for node ID {}'.format(_id))

			added_to_queue = 0
			parent_chem_doc = tree_dict[_id] # copy to overwrite later
			parent_chem_prod_of = parent_chem_doc['prod_of']

			# Assign unique number
			for (rxn, mols) in children:
				rxn_id = current_id.value
				current_id.value += 1 # this is only okay because there is ONE coordinator process

				# For the parent molecule, record child reactions
				parent_chem_prod_of.append(rxn_id)
			
				# For the reaction, keep track of children IDs
				chem_ids = []
				for mol in mols:
					
					# New chemical?
					if mol not in chem_to_id:

						chem_id = current_id.value
						current_id.value += 1 # this is only okay because there is ONE coordinator process

						# Check if buyable
						ppg = pricer.lookup_smiles(mol, alreadyCanonical = True)

						tree_dict[chem_id] = {
							'smiles': mol,
							'prod_of': [],
							'rct_of': [rxn_id],
							'depth': parent_chem_doc['depth'] + 1,
							'ppg': ppg
						}
						chem_to_id[mol] = chem_id

						if ppg:
							#print('{} buyable!'.format(mol))
							buyable_leaves.append(chem_id)

						else:
							# Add to queue to get expanded
							try:
								expansion_queues[parent_chem_doc['depth']].put((chem_id, mol))
								added_to_queue += 1
								#print('Put {} (ID {}) in expansion_queue {}'.format(mol, chem_id, parent_chem_doc['depth']))
							except IndexError:
								# Maximum depth reached
								#print('Reached maximum depth, so will not expand around {}'.format(tree_dict[chem_id]))
								pass

					else:
						chem_id = chem_to_id[mol]

						# Overwrite this chemical node to record it is a reactant of this rxn
						chem_doc = tree_dict[chem_id]
						chem_doc['rct_of'] += [rxn_id]
						tree_dict[chem_id] = chem_doc

					# Save ID
					chem_ids.append(chem_id)

				# Record by overwriting the whole dict value
				rxn['rcts'] = chem_ids
				rxn['prod'] = _id
				rxn['depth'] = parent_chem_doc['depth'] + 0.5
				tree_dict[rxn_id] = rxn

			# Overwrite dictionary entry for the parent
			parent_chem_doc['prod_of'] = parent_chem_prod_of
			tree_dict[_id] = parent_chem_doc

			# Report
			#print('Added {} children to expansion queue from parent node ID {}'.format(added_to_queue, _id))

		except VanillaQueue.Empty:
			pass


class TreeBuilder:
	'''
	The Transformer class defines an object which can be used to perform
	one-step retrosyntheses for a given molecule.

	All of the information ends up in self.tree_dict, which contains unique keys for
	each chemical and each reaction. Those have information about their children as well.
	'''

	def __init__(self, nb_workers = 2, max_depth = 3, Pricer = None, RetroTransformer = None):

		# Number of workers to expand nodes
		self.nb_workers = nb_workers
		self.max_depth = max_depth
		self.target = None
		self.isRunning = False

		self.Pricer = Pricer
		self.RetroTransformer = RetroTransformer



	def start_building(self, smiles):
		'''
		Begin building the network
		'''

		self.isRunning = True
		self.target = smiles 

		# Expansion queues of nodes that have yet to be expanded - one for each depth
		# This allows us to do the equivalent of a DFS for expansion
		self.expansion_queues = [Queue() for i in range(self.max_depth - 1)]

		# Queue of results that have yet to be processed and IDed
		self.results_queue = Queue()

		# Dictionary storing the overall tree
		self.manager = Manager()
		self.tree_dict = self.manager.dict()
		self.chem_to_id = self.manager.dict()
		self.paused = self.manager.Value('i', 0)
		self.done = self.manager.Value('i', 0)
		self.current_id = self.manager.Value('i', 1)
		self.buyable_leaves = self.manager.list() # keep track of buyable leaves for searching later
		self.workers = []

		# Initialize the highest-priority queue with the target
		# it's okay that this is really a depth-0 chemical, since it will be the expansion grabbed
		self.current_id.value = 2
		self.expansion_queues[-1].put((1, smiles))
		self.tree_dict[1] = {
			'smiles': smiles,
			'prod_of': [],
			'rct_of': [],
			'depth': 0,
		}

		# Begin processes
		for i in range(self.nb_workers):
			p = Process(target = expansion_worker, args = (i, self.expansion_queues, self.results_queue, self.paused, self.done, self.RetroTransformer))
			self.workers.append(p)
			p.start()

		self.coordinator = Process(target = coordinator, args = (self.current_id, self.tree_dict, self.expansion_queues, self.results_queue, self.chem_to_id, self.paused, self.done, self.buyable_leaves, self.Pricer))
		self.coordinator.start()


	def stop_building(self, timeout = 15):
		'''
		Stop building
		'''
		
		# Signal done
		self.done.value = 1

		print('Changed `done` signal to True')
		print('giving workers {} seconds to terminate'.format(timeout))

		# Join up workers
		time.sleep(timeout)	
		for p in self.workers + [self.coordinator]:
			if p.is_alive():
				print('Process is still alive?')
				p.terminate()

		print('Stopped building, all processes done')
		self.isRunning = False


	def is_paused(self):
		return (self.pause.value == 1)

	def is_running(self):
		return self.isRunning

	def pause(self):
		self.paused.value = 1

	def is_target(self, smiles):
		return self.target == smiles

	def unpause(self):
		self.paused.value = 0

	def info_string(self):
		num_chemicals = 0
		num_reactions = 0
		at_depth = {}
		for node in self.tree_dict.values():
			depth = node['depth']
			if depth % 1 == 0:
				num_chemicals += 1
			else:
				num_reactions += 1
			if depth not in at_depth:
				at_depth[depth] = 1
			else:
				at_depth[depth] += 1
		string = 'The current network contains {} nodes ({} chemicals and {} reactions)<br/>\n'.format(num_chemicals + num_reactions, num_chemicals, num_reactions)
		for (depth, count) in sorted(at_depth.iteritems(), key = lambda x: x[0]):
			label = 'chemicals' if depth % 1 == 0 else 'reactions'
			string += '...at depth {:.1f}, {} {}<br/>\n'.format(depth, count, label)
		return string

	def get_trees(self):
		'''
		This function recursively generates plausible synthesis trees where each chemical is
		buyable. It is build as a generator, so it should be called multiple times until it returns
		None or an error
		'''

		def chem_dict(_id, smiles, ppg, children = []):
			return {
				'id': _id,
				'is_chemical': True,
				'smiles' : smiles,
				'ppg' : ppg,
				'children': children,
			}

		def rxn_dict(_id, info, children = []):
			return {
				'id': _id,
				'is_reaction': True,
				'info': info,
				'children': children,
			}

		return [chem_dict(1, self.target, 0, children = [
					rxn_dict(2, 'rxn1', children = [
						chem_dict(3, 'CCCO', 1),
						chem_dict(4, 'CCC[Br]', 1)
					])
				])]

		# Recursively generate paths
		for path in viable_paths_starting_from_chemical(source_id = 1):
			return chem_dict(source_id, self.tree_dict[source_id]['smiles'], children = path)

		
if __name__ == '__main__':

	from pymongo import MongoClient
	db_client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	TRANSFORM_DB = db_client['reaxys']['transforms_retro_v4']
	CHEMICAL_DB = db_client['reaxys']['chemicals']
	BUYABLE_DB = db_client['reaxys']['buyables']

	import makeit.retro.transformer as transformer 
	RetroTransformer = transformer.Transformer()
	mincount_retro = 500
	print('Loading retro templates...')
	RetroTransformer.load(TRANSFORM_DB, mincount = mincount_retro, get_retro = True, get_synth = False)
	print('Loaded {} retro templates'.format(RetroTransformer.num_templates))


	print('Loading prices...')
	import makeit.retro.pricer as pricer
	Pricer = pricer.Pricer()
	Pricer.load(CHEMICAL_DB, BUYABLE_DB)
	print('Loaded known prices')

	from rdkit import RDLogger
	lg = RDLogger.logger()
	lg.setLevel(RDLogger.CRITICAL)

	builder = TreeBuilder(Pricer = Pricer, RetroTransformer = RetroTransformer)
	builder.start_building('CCCCOCCC')
	print('Will wait 30 seconds')
	time.sleep(5)
	builder.stop_building(timeout = 10)

	print(builder.info_string())

	print(builder.get_trees())