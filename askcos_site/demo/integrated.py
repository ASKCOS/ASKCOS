# Atropine
TARGET = 'CN3[C@H]1CC[C@@H]3C[C@@H](C1)OC(=O)C(CO)c2ccccc2'
# Fluconazole
TARGET = 'OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F'

expansion_time = 45
max_depth = 5
max_branching = 5
max_trees = 25


import time
from pymongo import MongoClient
db_client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)

RETRO_TRANSFORMS = {
    'database': 'reaxys',
    'collection': 'transforms_retro_v4', # 'lowe' or 'chematica'
    'mincount': 500,
    'parallel': False,
    'nb_workers': 1,
}
SYNTH_TRANSFORMS = {
    'database': 'reaxys',
    'collection': 'transforms_forward_v1'   ,
    'mincount': 100, 
}
INSTANCES = {
    'database': 'reaxys',
    'collection': 'instances',
}
REACTIONS = {
    'database': 'reaxys',
    'collection': 'reactions',
}
CHEMICALS = {
    'database': 'reaxys',
    'collection': 'chemicals',
}
BUYABLES = {
    'database': 'reaxys',
    'collection': 'buyables',
}
SOLVENTS = {
    'database': 'reaxys',
    'collection': 'solvents',
}

PREDICTOR = {
    'trained_model_path': '/home/ccoley/ML Chemistry/Make-It/makeit/predict/output/reloading',
    'info': '01-23-17, model trained on 80k Reaxys examples, validated on 10k, tested on 10k. Nh1_200, Nh2_200, Nh3_200, l2_0, Nc_5000, enh_weight_0d1, context_weight_50, opt_adadelta, batch_5, moreFeatures'
}


database = db_client[RETRO_TRANSFORMS['database']]
RETRO_DB = database[RETRO_TRANSFORMS['collection']]
import makeit.retro.transformer as transformer 
RetroTransformer = transformer.Transformer(
	parallel = RETRO_TRANSFORMS['parallel'], 
	nb_workers = RETRO_TRANSFORMS['nb_workers'],
)
print('Initialized retrotransformer')
mincount_retro = RETRO_TRANSFORMS['mincount']
RetroTransformer.load(RETRO_DB, mincount = mincount_retro, get_retro = True, get_synth = False)
print('Loaded {} retro templates'.format(RetroTransformer.num_templates))
RETRO_FOOTNOTE = 'Using {} retrosynthesis templates (mincount {}) from {}/{}'.format(RetroTransformer.num_templates,
	mincount_retro, RETRO_TRANSFORMS['database'], RETRO_TRANSFORMS['collection'])

### Forward transformer 
database = db_client[SYNTH_TRANSFORMS['database']]
SYNTH_DB = database[SYNTH_TRANSFORMS['collection']]
SynthTransformer = transformer.Transformer()
mincount_synth = SYNTH_TRANSFORMS['mincount']
SynthTransformer.load(SYNTH_DB, mincount = 100000000000000, get_retro = False, get_synth = True)
print('Loaded {} forward templates'.format(SynthTransformer.num_templates))
SYNTH_FOOTNOTE = 'Using {} forward templates (mincount {}) from {}/{}'.format(SynthTransformer.num_templates,
	mincount_synth, SYNTH_TRANSFORMS['database'], SYNTH_TRANSFORMS['collection'])

### Databases
db = db_client[REACTIONS['database']]
REACTION_DB = db[REACTIONS['collection']]
RETRO_LIT_FOOTNOTE = 'Searched {} known reactions from literature'.format(REACTION_DB.count())

db = db_client[INSTANCES['database']]
INSTANCE_DB = db[INSTANCES['collection']]

db = db_client[CHEMICALS['database']]
CHEMICAL_DB = db[CHEMICALS['collection']]

db = db_client[BUYABLES['database']]
BUYABLE_DB = db[BUYABLES['collection']]

db = db_client[SOLVENTS['database']]
SOLVENT_DB = db[SOLVENTS['collection']]

### Prices
print('Loading prices...')
import makeit.retro.pricer as pricer
Pricer = pricer.Pricer()
Pricer.load(CHEMICAL_DB, BUYABLE_DB)
print('Loaded known prices')

# Builder
from makeit.webapp.treeBuilder import TreeBuilder 
builder = TreeBuilder(Pricer = Pricer, RetroTransformer = RetroTransformer)

# Intelligent predictor
from makeit.webapp.forwardPredictor import ForwardPredictor 
predictor = ForwardPredictor(nb_workers = 2, TRANSFORM_DB = SYNTH_DB, SOLVENT_DB = SOLVENT_DB)
predictor.load_templates(mincount = mincount_synth)
predictor.load_model(PREDICTOR['trained_model_path'])
PREDICTOR_FOOTNOTE = 'Results generated using {} forward synthetic templates (mincount {}) from {}/{}, scored by a trained machine learning model: '.format(predictor.num_templates,	mincount_synth, SYNTH_TRANSFORMS['database'], SYNTH_TRANSFORMS['collection']) + PREDICTOR['info']

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

## Context recommendation
from askcos_site.functions.nnPredictor import nn_condition_predictor, lshf_nn, rxd_ids
NN = nn_condition_predictor()
NN.outputString = False

import rdkit.Chem as Chem

TARGET = Chem.MolToSmiles(Chem.MolFromSmiles(TARGET), isomericSmiles = False)
print('Target: {}'.format(TARGET))
builder.start_building(TARGET, max_depth = max_depth, max_branching = max_branching)
print('Will wait {} seconds'.format(expansion_time))
time.sleep(expansion_time)
print(builder.info_string())

builder.stop_building(timeout = 5)
print('###########################################################')
time.sleep(2)
print('IDDDFS trees:')
trees = builder.get_trees_iddfs(max_trees = max_trees)
print(trees)
print('{} trees'.format(len(trees)))


# Get context / evaluate?
context_dict = {} 
plausibility_dict = {}

def check_this_reaction(tree, this_is_the_target = False):
	'''Evaluates this reaction and all its children'''
	global context_dict 
	global plausibility_dict 

	if not tree['children']:
		print('No children to expand!')
		return 1.0

	products = [tree['smiles']]
	rxn = tree['children'][0]
	necessary_reagent = rxn['necessary_reagent']
	# TODO: Add reagent for necessary reagent! Hopefully this will be taken care of by NN
	# context, but not gauranteed. The necessary_reagent should be part of the context
	# recommendation system
	reactants = [child['smiles'] for child in rxn['children']]

	rxn_smiles = '.'.join(sorted(reactants)) + '>>' + products[0]

	if rxn_smiles in plausibility_dict:
		plausible = plausibility_dict[rxn_smiles]
	
	else:

		

		[T1, t1, y1, slvt1, rgt1, cat1] = NN.step_condition([reactants, products])
		if slvt1 and slvt1[-1] == '.': 
			slvt1 = slvt1[:-1]
		if rgt1 and rgt1[-1] == '.':
			rgt1 = rgt1[:-1]
		if cat1 and cat1[-1] == '.':
			cat1 = cat1[:-1]
		
		context_dict[rxn_smiles] = [T1, t1, y1, slvt1, rgt1, cat1]

		# Merge cat and reagent
		if rgt1 and cat1: 
			rgt1 = rgt1 + '.' + cat1 
		elif cat1:
			rgt1 = cat1

		if '.' in slvt1:
			slvt1 = slvt1.split('.')[0]


		error = predictor.set_context(T = T1, reagents = rgt1, solvent = slvt1)
		if error is not None:
			print('Recommended context: {}'.format([T1, t1, y1, slvt1, rgt1, cat1]))
			raise ValueError('NN recommended unrecognizable context!')


		# OPTION 2 - waits until completed
		predictor.run_in_foreground(reactants = '.'.join(reactants), intended_product = products[0], quit_if_unplausible = True)
		plausible = predictor.is_plausible()
		outcomes = predictor.return_top(n = 3)

		print('REACTANTS: {}'.format(reactants))
		print('PRODUCTS: {}'.format(products))
		print('necessary reagents: {}'.format(necessary_reagent))
		print('Recommended context: {}'.format([T1, t1, y1, slvt1, rgt1, cat1]))
		print('Plausible? {}'.format(plausible))
		print(outcomes)

		if plausible: 
			plausibility_dict[rxn_smiles] = float(predictor.products[products[0]])
		else:
			plausibility_dict[rxn_smiles] = 0.0

	# Look at children?
	if plausible: 
		all_plausible = True
		for child in rxn['children']:
			child_plausible = check_this_reaction(child)
			if not child_plausible:
				all_plausible = False 

		# Now report
		if all_plausible and this_is_the_target:
			print('################################################')
			print('Found completely plausible tree!')
			print(tree)
			raw_input('...pause...')
		elif this_is_the_target:
			print('Found tree with unfeasible children...thats not helpful')
	elif this_is_the_target:
		print('Found tree with unfeasible first step...thats not helpful')

	return plausible


for tree in trees:
	check_this_reaction(tree) # build up the dictionaries without using this_is_the_target

# Let all the children processes finish
for i in range(15):
	print('sleeping...')
	time.sleep(1)

print('Anything completely plausible?')
trees = builder.get_trees_iddfs(max_trees = max_trees)
for tree in trees:
	check_this_reaction(tree, this_is_the_target = True)

quit(1)