from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import rdkit.Chem as Chem          
from rdkit.Chem import AllChem
import numpy as np
import Queue as VanillaQueue
from multiprocessing import Pool, cpu_count, Process, Manager, Queue
from functools import partial # used for passing args to multiprocessing
import time
import os
import sys

from makeit.embedding.descriptors import edits_to_vectors
from makeit.predict.summarize_reaction_outcome import summarize_reaction_outcome


def softmax(x):
	"""Compute softmax values for each sets of scores in x."""
	e_x = np.exp(x - np.max(x))
	return e_x / e_x.sum()

def template_worker(i, workers_done, apply_queue, results_queue, done, F_atom, F_bond, atom_desc_dict):
	'''
	The target for a worker Process to apple one template to a molecule
	'''

	while True:
		# If done, stop
		if done.value:
			print('Worker {} saw done signal, terminating'.format(i))
			break

		# Grab something off the queue
		try:
			(template_rxn, reactants) = apply_queue.get(timeout = 0.1) # short timeout
			outcomes = template_rxn.RunReactants([reactants])
			if not outcomes: continue # no match

			candidate_list = []
			for j, outcome in enumerate(outcomes):
				outcome = outcome[0] # all products represented as single mol by transforms
				try:
					outcome.UpdatePropertyCache()
					Chem.SanitizeMol(outcome)
					[a.SetProp('molAtomMapNumber', a.GetProp('old_molAtomMapNumber')) \
						for a in outcome.GetAtoms() \
						if 'old_molAtomMapNumber' in a.GetPropsAsDict()]
						
				
					# Reduce to largest (longest) product only
					candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY)
					candidate_smiles = max(candidate_smiles.split('.'), key = len)
					outcome = Chem.MolFromSmiles(candidate_smiles)
						
					# Find what edits were made
					edits = summarize_reaction_outcome(reactants, outcome)

					# Remove mapping before matching
					[x.ClearProp('molAtomMapNumber') for x in outcome.GetAtoms() if x.HasProp('molAtomMapNumber')] # remove atom mapping from outcome

					# Overwrite candidate_smiles without atom mapping numbers
					candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY)

					# Add to ongoing list
					if (candidate_smiles, edits) not in candidate_list:
						candidate_list.append((candidate_smiles, edits))

				except Exception as e:
					continue	

			# Convert to tensors
			for (candidate_smiles, edits) in candidate_list:
				edit_h_lost_vec, edit_h_gain_vec, \
					edit_bond_lost_vec, edit_bond_gain_vec = edits_to_vectors(edits, reactants, atom_desc_dict = atom_desc_dict)

				x_h_lost = np.zeros((1, 1, len(edit_h_lost_vec), F_atom))
				x_h_gain = np.zeros((1, 1, len(edit_h_gain_vec), F_atom))
				x_bond_lost = np.zeros((1, 1, len(edit_bond_lost_vec), F_bond))
				x_bond_gain = np.zeros((1, 1, len(edit_bond_gain_vec), F_bond))

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

				x = [x_h_lost, x_h_gain, x_bond_lost, x_bond_gain]
				results_queue.put((candidate_smiles, x))
				#print('Worker {} added candidate {} to results for scoring'.format(i, candidate_smiles))

				sys.stdout.flush()
				sys.stderr.flush()
			
		except VanillaQueue.Empty:
			workers_done[i] = True 
			#print('Worker found empty queue')
			time.sleep(1)
			pass


def coordinator(workers_done, coordinator_done, apply_queue, results_queue, done, model, contexts, products, intended_product, plausible, quit_if_unplausible):
	'''
	The coordinator
	'''

	intended_found = [False for i in range(len(contexts))]
	intended_score = [-9999999 for i in range(len(contexts))]
	max_score = [-9999999 for i in range(len(contexts))]
	max_product = ['' for i in range(len(contexts))]

	while True:
		# If done, stop
		if done.value:
			print('Coordinator saw done signal, stopping')
			break
		try:
			(product_smiles, x) = results_queue.get(timeout = 0.1)

			for i, context in enumerate(contexts):
				score = model.predict(x + context)[0][0]
				#print('Score for {}, {}'.format(product_smiles, score))
				
				# Save score in shared dictionary - score is a list though
				# must overwrite WHOLE list as shared dict value
				if products.has_key(product_smiles):
					prev_scores = products[product_smiles]
					prev_scores[i] = max(prev_scores[i], score)
					products[product_smiles] = prev_scores
				else:
					prev_scores = [-99999 for j in range(len(contexts))]
					prev_scores[i] = score
					products[product_smiles] = prev_scores

				# Save max
				if score > max_score[i]:
					max_score[i] = score
					max_product[i] = product_smiles
					#print('Current most-likely product for context {}: {}'.format(i, max_product[i]))

				# Is this the intended product?
				if product_smiles == intended_product:
					intended_found[i] = True 
					intended_score[i] = products[product_smiles][i]

		except VanillaQueue.Empty:
			pass

		# Is the intended product not the major product?
		definitely_unplausible = [] 
		for i in range(len(contexts)):
			definitely_unplausible.append(
				intended_found[i] and intended_score[i] + 0.693 < max_score[i]
			)
		if all(definitely_unplausible):
			plausible.value = 0 # should be zero by default, but just to be safe
			if quit_if_unplausible:
				print('### intended product {} ({}) is significantly less likely than the current maximum-likelihood product {} ({})'.format(intended_product, intended_score, max_product, max_score))
				print('Quitting and flushing apply_queue')
				# Flush queue
				while not apply_queue.empty(): apply_queue.get(timeout = 0.1)
				coordinator_done.value = 1
				# Wait for workers to finish
				while not all(workers_done): time.sleep(0.1)
				# Flush queues again
				while not apply_queue.empty(): apply_queue.get(timeout = 0.1)
				while not results_queue.empty(): results_queue.get(timeout = 0.1)
				break

		# Are we finished?
		# note: this can be changed to stop early
		if all(workers_done) and results_queue.empty():
			print('All workers done and results queue is empty')
			print('Ready to quit')
			# # Flush
			# try:
			# 	while True:
			# 		apply_queue.get()
			# except VanillaQueue.Empty:
			# 	pass
			coordinator_done.value = 1
			definitely_plausible = [] 
			for i in range(len(contexts)):
				definitely_plausible.append(
					intended_found[i] and intended_score[i] + 0.693 > max_score[i]
				)
			if any(definitely_plausible):
				plausible.value = 1
			break


class ForwardPredictor:
	'''
	A class where a proposed reaction step can be evaluated for its likelihood of success
	'''

	def __init__(self, nb_workers = 10, TRANSFORM_DB = None, SOLVENT_DB = None):

		# Number of workers to expand nodes
		self.nb_workers = nb_workers

		# Status
		self.isRunning = False

		# Forwarad templates
		self.templates = []
		self.template_source = None

		# Queue of templates that need to be applied
		self.apply_queue = Queue()

		# Queue of results that have yet to be processed 
		self.results_queue = Queue()

		# List of workers
		self.workers = []

		# Done signal
		self.manager = Manager()
		self.done = self.manager.Value('i', 0)
		self.products = self.manager.dict()
		self.atom_desc_dict = self.manager.dict()
		self.intended_product = ''
		self.workers_done = self.manager.list([False for i in range(self.nb_workers)])
		self.coordinator_done = self.manager.Value('i', 0)
		self.plausible = self.manager.Value('i', 0)

		# DBs
		self.TRANSFORM_DB = TRANSFORM_DB
		self.SOLVENT_DB = SOLVENT_DB

		# Feature sizes
		mol = Chem.MolFromSmiles('[C:1][C:2]')
		(a, _, b, _) = edits_to_vectors((['1'],[],[('1','2',1.0)],[]), mol)
		self.F_atom = len(a[0])
		self.F_bond = len(b[0])

		# Begin processes
		for i in range(self.nb_workers):
			p = Process(target = template_worker, args = (i, self.workers_done, self.apply_queue, self.results_queue, self.done, self.F_atom, self.F_bond, self.atom_desc_dict))
			self.workers.append(p)
			p.start()


	def load_templates(self, mincount = 25):
		'''Load templates from a collection'''

		for document in self.TRANSFORM_DB.find({'count': {'$gt': mincount}}, ['_id', 'count', 'reaction_smarts']):
			# Skip if no reaction SMARTS
			if 'reaction_smarts' not in document: continue
			reaction_smarts = str(document['reaction_smarts'])
			if not reaction_smarts: continue

			# Define dictionary
			template = {
				'reaction_smarts': 		reaction_smarts,
				'count': 				document['count'] if 'count' in document else 0,
				'_id':	 				document['_id'] if '_id' in document else -1,
			}

			try:
				rxn_f = AllChem.ReactionFromSmarts('(' + reaction_smarts.replace('>>', ')>>(') + ')')
				#if rxn_f.Validate() == (0, 0):
				if rxn_f.Validate()[1] == 0:
					template['rxn_f'] = rxn_f
				else:
					continue
			except Exception as e:
				print('Couldnt load forward: {}: {}'.format(reaction_smarts, e))
				continue

			# Add to list
			self.templates.append(template)

		self.num_templates = len(self.templates)
		self.templates[:] = [x for x in sorted(self.templates, key = lambda z: z['count'], reverse = True)]

		print('Loaded {} templates'.format(self.num_templates))

	def load_model(self, folder):
		'''Load a model'''

		# Get model args
		ARGS_FPATH = os.path.join(folder, 'args.json')
		with open(ARGS_FPATH, 'r') as fid:
			import json
			args = json.load(fid)

		N_h2 = int(args['Nh2'])
		N_h1 = int(args['Nh1'])
		N_h3 = int(args['Nh3'])
		N_hf = int(args['Nhf'])
		l2v = float(args['l2'])
		lr = float(args['lr'])
		context_weight = float(args['context_weight'])
		enhancement_weight = float(args['enhancement_weight'])
		optimizer          = args['optimizer']
		inner_act          = args['inner_act']
		TARGET_YIELD       = False

		from makeit.predict.reaxys_score_candidates_from_edits_compact_v2 import build
		self.model = build(F_atom = self.F_atom, F_bond = self.F_bond, N_h1 = N_h1, 
				N_h2 = N_h2, N_h3 = N_h3, N_hf = N_hf, l2v = l2v, inner_act = inner_act,
				context_weight = context_weight, enhancement_weight = enhancement_weight, TARGET_YIELD = TARGET_YIELD,
				absolute_score = True)

		WEIGHTS_FPATH = os.path.join(folder, 'weights.h5')
		self.model.load_weights(WEIGHTS_FPATH, by_name = True)

	def stop(self, timeout = 15):
		'''
		Stop building
		'''

		if not self.is_running(): return
		
		# Signal done
		if not self.done.value:
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

	def set_context(self, reagents = '', solvent = '', T = ''):

		# Temperature is easy
		try:
			T = float(T)
		except TypeError:
			return 'Cannot convert temperature {} to float'.format(T)

		# Solvent needs a lookup
		solvent_mol = Chem.MolFromSmiles(solvent)
		if solvent_mol: 
			doc = self.SOLVENT_DB.find_one({'_id': Chem.MolToSmiles(solvent_mol)})
		else:
			doc = self.SOLVENT_DB.find_one({'name': solvent})
		if not doc:
			return 'Could not parse solvent {}'.format(solvent)
		solvent_vec = [doc['c'], doc['e'], doc['s'], doc['a'], doc['b'], doc['v']]
		solvent = doc['name']

		# Unreacting reagents
		reagents = [Chem.MolFromSmiles(reagent) for reagent in reagents.split('.')]
		if None in reagents:
			return 'Could not parse all reagents!'
		reagent_fp = np.zeros((1, 256))
		for reagent in reagents:
			reagent_fp += np.array(AllChem.GetMorganFingerprintAsBitVect(reagent, 2, nBits = 256))
		
		# Save list
		self.context = [[reagent_fp, np.reshape(np.array(solvent_vec), (1, 6)), np.reshape(np.array(T), (1, 1))]]

	def set_contexts(self, contexts):
		'''Set multiple contexts to try at once'''

		self.contexts = []
		self.context_labels = []
		errors = []
		for (T, reagents, solvent) in contexts:
			# Temperature is easy
			try:
				T = float(T)
			except TypeError:
				#return 'Cannot convert temperature {} to float'.format(T)
				errors.append('cannot convert temp {} to float'.format(T))
				continue

			# Solvent needs a lookup
			solvent_mol = Chem.MolFromSmiles(solvent)
			if solvent_mol: 
				doc = self.SOLVENT_DB.find_one({'_id': Chem.MolToSmiles(solvent_mol)})
			else:
				doc = self.SOLVENT_DB.find_one({'name': solvent})
			if not doc:
				#return 'Could not parse solvent {}'.format(solvent)
				errors.append('could not parse solvent {}'.format(solvent))
				continue
			solvent_vec = [doc['c'], doc['e'], doc['s'], doc['a'], doc['b'], doc['v']]
			solvent = doc['name']

			# Unreacting reagents
			reagents_mols = [Chem.MolFromSmiles(reagent) for reagent in reagents.split('.')]
			if None in reagents_mols:
				#return 'Could not parse all reagents!'
				errors.append('could not parse all reagents')
				continue
			reagent_fp = np.zeros((1, 256))
			for reagent in reagents_mols:
				reagent_fp += np.array(AllChem.GetMorganFingerprintAsBitVect(reagent, 2, nBits = 256))
			
			# Save list
			self.contexts.append([reagent_fp, np.reshape(np.array(solvent_vec), (1, 6)), np.reshape(np.array(T), (1, 1))])
			self.context_labels.append('T:{}, rgt:{}, solv:{}'.format(T, reagents, solvent))
			errors.append(False)
		return errors

	def run(self, reactants = '', intended_product = '', quit_if_unplausible = False, mincount = 0):
		
		# Force clear products
		del self.products 
		self.products = self.manager.dict()

		reactants = Chem.MolFromSmiles(reactants)
		if not reactants: 
			print('Could not parse reactants {}'.format(reactants))
			raise ValueError('Could not parse reactants')
		print('Number of reactant atoms: {}'.format(len(reactants.GetAtoms())))
		# Report current reactant SMILES string
		[a.ClearProp('molAtomMapNumber') for a in reactants.GetAtoms() if a.HasProp('molAtomMapNumber')]
		print('Reactants w/o map: {}'.format(Chem.MolToSmiles(reactants)))
		# Add new atom map numbers
		[a.SetProp('molAtomMapNumber', str(i+1)) for (i, a) in enumerate(reactants.GetAtoms())]
		# Report new reactant SMILES string
		print('Reactants w/ map: {}'.format(Chem.MolToSmiles(reactants)))

		self.reactants_smiles = Chem.MolToSmiles(reactants)
		self.intended_product = intended_product
		self.quit_if_unplausible = quit_if_unplausible

		# Pre-calc descriptors for this set of reactants
		for (key, val) in edits_to_vectors([], reactants, return_atom_desc_dict = True).iteritems():
			self.atom_desc_dict[key] = val

		# Clear queues (just to be safe)
		while not self.apply_queue.empty(): self.apply_queue.get()
		while not self.results_queue.empty(): self.results_queue.get()

		# Load up queue
		for template in self.templates:
			if template['count'] < mincount: break # dont queue up
			self.apply_queue.put((template['rxn_f'], reactants))
		for i in range(self.nb_workers):
			self.workers_done[i] = False

		# Start the coordinator
		self.coordinator_done.value = 0
		self.plausible.value = 0
		self.coordinator = Process(target = coordinator, args = (self.workers_done, self.coordinator_done, self.apply_queue, self.results_queue, self.done, self.model, self.contexts, self.products, self.intended_product, self.plausible, self.quit_if_unplausible))
		self.coordinator.start()
		self.isRunning = True

		print('Added forward templates to queue and started coordinator')

	def is_running(self):
		return bool(self.isRunning)

	def is_coord_done(self):
		return bool(self.coordinator_done.value)

	def is_done(self):
		return bool(self.done.value)

	def is_plausible(self):
		return bool(self.plausible.value)

	def run_in_foreground(self, **kwargs):
		self.run(**kwargs)
		while not self.is_coord_done(): pass

	def num_prods(self):
		return self.products.__len__()

	def return_top(self, n = 25):
		'''Return the top n outcomes'''

		all_outcomes = [] 
		for j, context in enumerate(self.contexts):

			# Pull the score for this specific context
			sorted_prods = sorted(
				[(prod, score[j]) for (prod, score) in self.products.items()], 
				key = lambda x: x[1], reverse = True
			)
			# print(sorted_prods)
			if not sorted_prods: return []

			prods, scores = zip(*sorted_prods)
			probs = softmax(scores)
			probs_truncated = softmax(scores[:n])

			outcomes = []
			for i in range(len(prods)):
				if i == n: break 
				if scores[i] < -500: break
				outcomes.append({
					'rank': i + 1,
					'smiles': prods[i],
					'score': '{:.2f}'.format(scores[i]),
					'prob': '{:.2e}'.format(probs[i]),
					'prob_trunc': '{:.2e}'.format(probs_truncated[i]),
				})

			all_outcomes.append(outcomes)
		return all_outcomes

if __name__ == '__main__':

	from pymongo import MongoClient
	db_client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	TRANSFORM_DB = db_client['reaxys']['transforms_forward_v1']
	SOLVENT_DB = db_client['reaxys']['solvents']

	from rdkit import RDLogger
	lg = RDLogger.logger()
	lg.setLevel(RDLogger.CRITICAL)

	predictor = ForwardPredictor(nb_workers = 2, TRANSFORM_DB = TRANSFORM_DB, SOLVENT_DB = SOLVENT_DB)
	predictor.load_templates(mincount = 100)
	predictor.load_model('/home/ccoley/ML Chemistry/Make-It/makeit/predict/output/reloading')
	predictor.set_contexts([
		('20', '[Na+].[O-]', 'water'),
		('40', '[K+].[O-]', 'water'),
		('100', '', 'toluene'),
	])

	# OPTION 1
	# predictor.run(reactants = 'CCCCCCCO.CC[Br]', intended_product = 'CCCCCCCOCC')
	# while not predictor.is_coord_done(): pass
	# predictor.stop(timeout = 0.5)

	# OPTION 2 - waits until completed
	predictor.run_in_foreground(reactants = 'CCCCCCCO.CC[Br]', intended_product = 'CCCCCCC=O', quit_if_unplausible = True)

	plausible = predictor.is_plausible()
	print('Plausible? {}'.format(plausible))
	outcomes = predictor.return_top(n = 3)
	for i, outcome in enumerate(outcomes):
		print('### FOR CONTEXT {}'.format(predictor.context_labels[i]))
		print(outcome)



	predictor.run_in_foreground(reactants = 'CCCCCCCO.CC[Br]', intended_product = 'CCCCCCC=O', quit_if_unplausible = False)

	plausible = predictor.is_plausible()
	print('Plausible? {}'.format(plausible))
	outcomes = predictor.return_top(n = 3)
	for i, outcome in enumerate(outcomes):
		print('### FOR CONTEXT {}'.format(predictor.context_labels[i]))
		print(outcome)

	predictor.stop(timeout = 0.5)