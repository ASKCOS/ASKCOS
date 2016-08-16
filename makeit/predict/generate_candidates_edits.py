# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import argparse
import numpy as np     	      	   # for simple calculations
import os                          # for saving
import sys
import rdkit.Chem as Chem
import makeit.retro.transformer as transformer 
from makeit.retro.canonicalization import SmilesFixer
from pymongo import MongoClient    # mongodb plugin
from makeit.predict.summarize_reaction_outcome import summarize_reaction_outcome
from makeit.embedding.descriptors import edits_to_vectors, oneHotVector # for testing
import re
import time
from tqdm import tqdm

def main(template_collection = 'lowe_refs_general_v3', reaction_collection = 'lowe_1976-2013_USPTOgrants_reactions',
		candidate_collection = 'candidate_edits_mincount50', USE_REACTIONSMILES = True, mincount = 4, n_max = 50, seed = None, 
		outfile = '.', singleonly = True, check = True, log = False):

	from rdkit import RDLogger
	lg = RDLogger.logger()
	lg.setLevel(4)

	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['askcos_transforms']
	templates = db[template_collection]
	Transformer = transformer.Transformer()
	Transformer.load(templates, mincount = mincount, get_retro = False, get_synth = True)
	print('Out of {} database templates,'.format(templates.count()))
	print('Loaded {} templates'.format(Transformer.num_templates))
	Transformer.reorder()
	print('Sorted by count, descending')

	db = client['reaction_examples']
	reactions = db[reaction_collection]

	db = client['prediction']
	candidates = db[candidate_collection]

	if check:
		done_ids = [doc['reaction_id'] for doc in candidates.find(
			{'reaction_collection': reaction_collection, 'found': True}
		)]
		print('Checked completed entries')

	else:
		done_ids = []

	# Define generator
	class Randomizer():
		def __init__(self, seed, done_ids = []):
			self.done_ids = done_ids
			np.random.seed(seed)
			if outfile:
				with open(os.path.join(outfile, 'seed.txt'), 'w') as fid:
					fid.write('{}'.format(seed))
		def get_rand(self):
			'''Random WITHOUT replacement'''
			while True:
				doc = reactions.find({'products': {'$size': 1}, \
					'products.0.smiles': {'$not': re.compile('.*\..*')},
					'random': { '$gte': np.random.random()}}).sort('random', 1).limit(1)[0]
				# Check if it was done before this script ran
				if doc['_id'] in self.done_ids: 
					continue
				# Check if it has been done since (when running multiple instances in parallel)
				if candidates.count({'reaction_id': doc['_id'], 
						'reaction_collection': reaction_collection}) > 0:
					continue
				self.done_ids.append(doc['_id'])
				yield doc

	if seed == None:
		seed = np.random.randint(100000)
	else:
		seed = int(seed)
	randomizer = Randomizer(seed, done_ids = done_ids)
	generator = enumerate(randomizer.get_rand())

	smilesfixer = SmilesFixer()

	# LOGGING
	if log:
		flog = open('GENERATE_CANDIDATES_LOG_{}.txt'.format(seed), 'w')
		flog.write('mincount: {}\n'.format(mincount))
		flog.write('number of templates: {}\n'.format(Transformer.num_templates))

	try:
		for i, reaction in generator:

			# LOGGING
			start_time = time.time()

			if i == n_max: 
				break

			print('#########')
			print('## RXN {}'.format(i))
			print('#########')

			if bool(USE_REACTIONSMILES):
				rxn_smiles = reaction['reaction_smiles'].split(' ')[0]
				all_smiles = [smilesfixer.fix_smiles(x) for x in rxn_smiles.split('>')[0].split('.')]
				mol = Chem.MolFromSmiles(rxn_smiles.split('>')[2])
			else:
				all_smiles =  [smilesfixer.fix_smiles(x['smiles']) for x in reaction['reactants']]
				if 'catalysts' in reaction:
					all_smiles += [smilesfixer.fix_smiles(x['smiles']) for x in reaction['catalysts']] 
				if 'spectators' in reaction:
					all_smiles += [smilesfixer.fix_smiles(x['smiles']) for x in reaction['spectators']] 
				mol = Chem.MolFromSmiles(reaction['products'][0]['smiles'])
			
			# Define target (true) product smiles
			[x.ClearProp('molAtomMapNumber') for x in mol.GetAtoms()] # remove atom mapping from target prod
			target_smiles = smilesfixer.fix_smiles(Chem.MolToSmiles(mol, isomericSmiles = USE_STEREOCHEMISTRY))
			if singleonly:
				target_smiles = max(target_smiles.split('.'), key = len)

			# Load reactant molecules
			reactants = Chem.MolFromSmiles('.'.join(all_smiles))
			n_reactant_atoms = len(reactants.GetAtoms())
			print('Number of reactant atoms: {}'.format(n_reactant_atoms))
			if n_reactant_atoms > 100:
				print('Skipping huge molecule! N_reactant_atoms = {}'.format(n_reactant_atoms))
				continue

			# Report current reactant SMILES string
			print('Reactants w/o map: {}'.format(Chem.MolToSmiles(reactants)))
			# Add new atom map numbers
			[a.SetProp('molAtomMapNumber', str(i+1)) for (i, a) in enumerate(reactants.GetAtoms())]
			# Report new reactant SMILES string
			print('Reactants w/ map: {}'.format(Chem.MolToSmiles(reactants)))

			found_true = False
			candidate_edits = [] # list of tuples of (product smiles, edits required)
			for template in tqdm(Transformer.templates):
				# Perform transformation
				try:
					outcomes = template['rxn_f'].RunReactants([reactants])
				except Exception as e:
					if v: print(e)
					continue
				if not outcomes: continue # no match
				for j, outcome in enumerate(outcomes):
					outcome = outcome[0] # all products represented as single mol by transforms
					try:
						outcome.UpdatePropertyCache()
						Chem.SanitizeMol(outcome)
						[a.SetProp('molAtomMapNumber', a.GetProp('old_molAtomMapNumber')) \
							for (i, a) in enumerate(outcome.GetAtoms()) \
							if 'old_molAtomMapNumber' in a.GetPropsAsDict()]
					except Exception as e:
						if v: print(e)
						continue
					if v: print('Outcome SMILES: {}'.format(Chem.MolToSmiles(outcome)))

					
					
					# Reduce to largest (longest) product only?
					candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY)
					if singleonly: 
						candidate_smiles = max(candidate_smiles.split('.'), key = len)
						outcome = Chem.MolFromSmiles(candidate_smiles)
						
					# Find what edits were made
					edits = summarize_reaction_outcome(reactants, outcome)
					if v: print(edits)

					# Remove mapping before matching
					[x.ClearProp('molAtomMapNumber') for x in outcome.GetAtoms() if x.HasProp('molAtomMapNumber')] # remove atom mapping from outcome
					if Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY) == target_smiles:
						if v: print('Matched true [{}]'.format(target_smiles))
						found_true = True

					# Overwrite candidate_smiles without atom mapping numbers
					candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles = USE_STEREOCHEMISTRY)

					# Add to ongoing list
					if (candidate_smiles, edits) not in candidate_edits:
						candidate_edits.append((candidate_smiles, edits))

			# Prepare doc and insert
			if found_true: print('Found true product')
			doc = {
				'_id': reaction['_id'],
				'reaction_collection': reaction_collection,
				'reaction_id': reaction['_id'],
				'reactant_smiles': Chem.MolToSmiles(reactants, isomericSmiles = USE_STEREOCHEMISTRY),
				'edit_candidates': candidate_edits,
				'product_smiles_true': target_smiles,
				'found': found_true,
				'num_candidates': len(candidate_edits),
			}
			try:
				res = candidates.insert(doc)
			except Exception as e:
				print(e)
				continue

			# LOGGING
			end_time = time.time()
			print('time: {}'.format(end_time - start_time))
			print('unique edit sets using longest prod: {}'.format(len(candidate_edits)))
			if log: flog.write('{}\t{}\t{}\t{}\n'.format(i, n_reactant_atoms, len(candidate_edits), end_time - start_time))

			#raw_input('Pause')

	except Exception as e:
		print('Error! {}'.format(e))

	if log: flog.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing; defaults to False')
	parser.add_argument('-n', '--num', type = int, default = 50,
						help = 'Maximum number of records to examine; defaults to 50')
	parser.add_argument('--reaction_collection', type = str, default = 'lowe_1976-2013_USPTOgrants_reactions',
						help = 'Collection of reaction_examples to use; defaults to lowe_1976-2013_USPTOgrants_reactions')
	parser.add_argument('--template_collection', type = str, default = 'lowe_refs_general_v3',
						help = 'Collection of templates to use; defaults to lowe_refs_general_v3')
	parser.add_argument('--candidate_collection', type = str, default = 'candidate_edits_8_9_16', 
						help = 'Collection of candidates to write to; defaults to candidate_edits')
	parser.add_argument('--seed', type = int, default = None,
						help = 'Seed for random number generator')
	parser.add_argument('--rxnsmiles', type = bool, default = True,
						help = 'Use reaction_smiles for species, not pre-parsed; defaults to true')
	parser.add_argument('--mincount', type = int, default = 50,
						help = 'Minimum template count to include in transforms; defaults to 50')
	parser.add_argument('--singleonly', type = bool, default = True,
						help = 'Whether to record major product only; defaults to True')
	parser.add_argument('--check', type = bool, default = True,
						help = 'Whether to check current collection to see if reaction example has been done')
	parser.add_argument('--log', type = bool, default = True,
						help = 'Whether to log wall times / number of candidate atoms / etc., default True')
	args = parser.parse_args()
	v = bool(args.v)

	main(template_collection = args.template_collection, reaction_collection = args.reaction_collection,
		candidate_collection = args.candidate_collection, seed = args.seed, n_max = int(args.num), 
		mincount = int(args.mincount), USE_REACTIONSMILES = bool(args.rxnsmiles), singleonly = bool(args.singleonly),
		check = bool(args.check), log = bool(args.log))
