import os
import sys
import time
import argparse
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

def main(TARGET, expansion_time=60.0, max_depth=5, max_branching=30, max_trees=1000):
    global context_dict
    global plausibility_dict
    global plausible_trees
    global done_trees
    global print_and_log

    TARGET = Chem.MolToSmiles(Chem.MolFromSmiles(TARGET), isomericSmiles=False)
    print_and_log('Target: {}'.format(TARGET), success=True)
    builder.start_building(TARGET, max_depth=max_depth, max_branching=max_branching)
    print_and_log('Will wait {} seconds'.format(expansion_time), success=True)
    time.sleep(expansion_time)
    print_and_log(builder.info_string(), success=True)

    builder.stop_building(timeout=5)
    print_and_log('###########################################################', success=True)
    time.sleep(2)

    # Function to recontruct the plausible tree
    def construct_plausible_tree(tree, depth=0):
        if not tree['children']:
            print_and_log('Depth {}: buyable chemical: {}'.format(depth, tree['smiles']), success=True)
        else:
            product = tree['smiles']
            rxn = tree['children'][0]
            necessary_reagent = rxn['necessary_reagent']
            # TODO: Add reagent for necessary reagent! Hopefully this will be taken care of by NN
            # context, but not gauranteed. The necessary_reagent should be part of the context
            # recommendation system
            reactants = [child['smiles'] for child in rxn['children']]

            rxn_smiles = '.'.join(sorted(reactants)) + '>>' + product
            print_and_log('Depth {}: reaction smiles: {} plausibility: {}'.format(
                depth, rxn_smiles, plausibility_dict[rxn_smiles]), success=True)
            print_and_log('Recommended context: {}'.format(context_dict[rxn_smiles]), success=True)
            print_and_log('Necessary reagents: {}'.format(necessary_reagent), success=True)
            depth += 1
            for child in rxn['children']:
                construct_plausible_tree(child, depth=depth)

    # Now define the reaction-checking function
    # Get context / evaluate?

    def check_this_reaction(tree, this_is_the_target=False):
        '''Evaluates this reaction and all its children'''
        print_and_log('---------------')

        if not tree['children']:
            print_and_log('No children to expand!')
            return 1.0

        products = [tree['smiles']]
        rxn = tree['children'][0]
        necessary_reagent = rxn['necessary_reagent']
        # TODO: Add reagent for necessary reagent! Hopefully this will be taken care of by NN
        # context, but not gauranteed. The necessary_reagent should be part of the context
        # recommendation system
        reactants = [child['smiles'] for child in rxn['children']]

        rxn_smiles = '.'.join(sorted(reactants)) + '>>' + products[0]
        print_and_log('reaction smiles: {}'.format(rxn_smiles))

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

            # Add reagent to reactants for the purpose of forward prediction
            # this is important for necessary reagents, e.g., chlorine source
            if rgt1:
                reactants.append(rgt1)

            # Merge cat and reagent
            if rgt1 and cat1:
                rgt1 = rgt1 + '.' + cat1
            elif cat1:
                rgt1 = cat1

            if '.' in slvt1:
                slvt1 = slvt1.split('.')[0]

            error = predictor.set_context(T=T1, reagents=rgt1, solvent=slvt1)
            if error is not None:
                print_and_log('Recommended context: {}'.format([T1, t1, y1, slvt1, rgt1, cat1]))
                print_and_log('NN recommended unrecognizable context for {}'.format(rxn_smiles))
                plausible = 0.

            else:
                # Run predictor
                predictor.run_in_foreground(reactants='.'.join(reactants), intended_product=products[0],
                                            quit_if_unplausible=True)
                outcomes = predictor.return_top(n=RANK_THRESHOLD_FOR_INCLUSION)
                plausible = 0.
                for outcome in outcomes:
                    if outcome['smiles'] == products[0]:
                        plausible = float(outcome['prob'])
                        break
                print_and_log('necessary reagents: {}'.format(necessary_reagent))
                print_and_log('Recommended context: {}'.format([T1, t1, y1, slvt1, rgt1, cat1]))
                print_and_log('plausible? {}'.format(plausible))
                print_and_log(outcomes[:min(3, len(outcomes))])

            plausibility_dict[rxn_smiles] = plausible

        # Look at children?
        if plausible:
            all_plausible = True
            for child in rxn['children']:
                child_plausible = check_this_reaction(child)
                if not child_plausible:
                    all_plausible = False

            # Now report
            if all_plausible and this_is_the_target:
                plausible_trees.append(tree)
                print_and_log('################################################', success=True)
                print_and_log('Found completely plausible tree!', success=True)
                print_and_log(tree, success=True)
                construct_plausible_tree(tree, depth=0)
                print_and_log('################################################', success=True)
            elif this_is_the_target:
                print_and_log('Found tree with unfeasible children...thats not helpful')
        elif this_is_the_target:
            print_and_log('Found tree with unfeasible first step...thats not helpful')

        return plausible


    # Define chemical and reaction helper functions (redundant with builder)
    def chem_dict(_id, smiles, ppg, children=[]):
        '''Chemical object as expected by website'''
        return {
            'id': _id,
            'is_chemical': True,
            'smiles': smiles,
            'ppg': ppg,
            'children': children,
        }


    def rxn_dict(_id, info, necessary_reagent='', num_examples=0, children=[]):
        '''Reaction object as expected by website'''
        return {
            'id': _id,
            'is_reaction': True,
            'info': info,
            'necessary_reagent': necessary_reagent,
            'num_examples': num_examples,
            'children': children,
        }


    # Now define the search (redundant with builder)
    def IDDFS_generator(tree_dict={}, max_depth=3):
        '''Get viable synthesis trees using an iterative deepening depth-first search'''

        def DLS_chem(chem_id, depth, headNode=False):
            '''Expand at a fixed depth for the current node chem_id.

            headNode indicates whether this is the first (head) node, in which case
            we should expand it even if it is itself buyable'''

            if depth == 0:
                # Not allowing deeper - is this buyable?
                if tree_dict[chem_id]['ppg']:
                    yield []  # viable node, calling function doesn't need children
            else:
                # Do we need to go deeper?
                if tree_dict[chem_id]['ppg'] and not headNode:
                    yield []  # Nope, this is a viable node
                else:
                    # Try going deeper via DLS_rxn function
                    for rxn_id in tree_dict[chem_id]['prod_of']:
                        rxn_info_string = ''
                        for path in DLS_rxn(rxn_id, depth):
                            yield [rxn_dict(rxn_id, rxn_info_string, tree_dict[rxn_id]['necessary_reagent'],
                                            tree_dict[rxn_id]['num_examples'], children=path)]

        def DLS_rxn(rxn_id, depth):
            '''Return children paths starting from a specific rxn_id

            Weird case handling based on 1, 2, or 3 reactants'''

            # Define generators for each reactant node - decrement the depth!
            generators = [DLS_chem(chem_id, depth - 1) for chem_id in tree_dict[rxn_id]['rcts']]

            # Only one reactant? easy!
            if len(generators) == 1:
                chem_id = tree_dict[rxn_id]['rcts'][0]
                for path in generators[0]:
                    yield [
                        chem_dict(chem_id, tree_dict[chem_id]['smiles'], tree_dict[chem_id]['ppg'], children=path)
                    ]

            # Two reactants? want to capture all combinations of each node's options
            elif len(generators) == 2:
                chem_id0 = tree_dict[rxn_id]['rcts'][0]
                chem_id1 = tree_dict[rxn_id]['rcts'][1]
                for path0 in generators[0]:
                    for path1 in generators[1]:
                        yield [
                            chem_dict(chem_id0, tree_dict[chem_id0]['smiles'], tree_dict[chem_id0]['ppg'],
                                      children=path0),
                            chem_dict(chem_id1, tree_dict[chem_id1]['smiles'], tree_dict[chem_id1]['ppg'],
                                      children=path1)
                        ]

            # Three reactants? This is not elegant...
            elif len(generators) == 3:
                chem_id0 = tree_dict[rxn_id]['rcts'][0]
                chem_id1 = tree_dict[rxn_id]['rcts'][1]
                chem_id2 = tree_dict[rxn_id]['rcts'][2]
                for path0 in generators[0]:
                    for path1 in generators[1]:
                        for path2 in generators[2]:
                            yield [
                                chem_dict(chem_id0, tree_dict[chem_id0]['smiles'], tree_dict[chem_id0]['ppg'],
                                          children=path0),
                                chem_dict(chem_id1, tree_dict[chem_id1]['smiles'], tree_dict[chem_id1]['ppg'],
                                          children=path1),
                                chem_dict(chem_id2, tree_dict[chem_id2]['smiles'], tree_dict[chem_id2]['ppg'],
                                          children=path2),
                            ]

            else:
                print_and_log('Too many reactants! Only have cases 1-3 programmed')
                raise ValueError('Too many reactants! Only have cases 1-3 programmed')

        for depth in range(max_depth):
            print_and_log('~~ looking at depth {}'.format(depth))
            try:
                for path in DLS_chem(1, depth, headNode=True):
                    yield chem_dict(1, tree_dict[1]['smiles'], tree_dict[1]['ppg'], children=path)
            except Exception as e:
                print_and_log('Path with exception: {}'.format(e))
                yield None


    # Look for trees
    print_and_log('Looking for trees now')
    import hashlib
    import json

    trees = []
    counter = 0
    for tree in IDDFS_generator(
            tree_dict=builder.tree_dict,
            max_depth=builder.max_depth):

        # Step 0 - check if with exception
        if tree is None:
            print_and_log('Found tree with exceptions, skipping')
            continue

        # Step 1 - check if unique
        hashkey = hashlib.sha1(json.dumps(tree, sort_keys=True)).hexdigest()
        if hashkey in done_trees:
            print_and_log('Found duplicate tree, skipping')
            continue

        # Step 2 - check feasibility
        check_this_reaction(tree, this_is_the_target=True)

        # Added into the done_trees set
        done_trees.add(hashkey)

        counter = counter + 1
        if counter == max_trees:
            print('Reach max number of trees: {}'.format(max_trees))
            break

    print_and_log('Got out of final IDDFS_generator loop, done with program')
    builder.stop_building()
    predictor.stop()
    print_and_log('Stopped builder and predictor')
    # del builder
    # del predictor
    # del NN
    fid.close()
    fid_success.close()

    # quit(1)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--TARGET', type=str,
                        help='Name of chemical')
    parser.add_argument('--expansion_time', type=float, default=60.0,
                        help='Retrosynthesis tree expansion time, default 60.0 sec')
    parser.add_argument('--max_depth', type=int, default=3,
                        help='Retrosynthesis tree expansion depth, default 3')
    parser.add_argument('--max_branching', type=int, default=30,
                        help='Retrosynthesis tree expansion branching ratio, default 30')
    parser.add_argument('--max_trees', type=int, default=1000,
                        help='Maxinum number of trees to examine, default 1000')
    parser.add_argument('--retro_mincount', type=int, default=100,
                        help = 'Minimum template count for retrosynthetic templates, default 100')
    parser.add_argument('--synth_mincount', type=int, default=50,
                        help = 'Minimum template count for synthetic templates, default 50')
    parser.add_argument('--rank_threshold', type=int, default=10,
                        help = 'Rank threshold for considering a reaction to be plausible, default 10')
    parser.add_argument('--max_expansions', type=int, default=8,
                        help = 'Maximum number of search expansions to try, default 8')
    args = parser.parse_args()

    

    expansion_time = float(args.expansion_time)
    max_depth = int(args.max_depth)
    max_branching = int(args.max_branching)
    max_trees = int(args.max_trees)
    TARGET = str(args.TARGET)
    RANK_THRESHOLD_FOR_INCLUSION = int(args.rank_threshold)

    import urllib2
    smiles = urllib2.urlopen('https://cactus.nci.nih.gov/chemical/structure/{}/smiles'.format(TARGET)).read()
    print('Resolved {} smiles -> {}'.format(TARGET, smiles))
    TARGET_LABEL = TARGET 
    TARGET = smiles

    if True: # lazy setup
        

        from pymongo import MongoClient

        db_client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)

        RETRO_TRANSFORMS = {
            'database': 'reaxys',
            'collection': 'transforms_retro_v4',  # 'lowe' or 'chematica'
            'mincount': int(args.retro_mincount),
            'parallel': False, # treebuilder ignores this
            'nb_workers': 1,
        }
        SYNTH_TRANSFORMS = {
            'database': 'reaxys',
            'collection': 'transforms_forward_v1',
            'mincount': int(args.synth_mincount),
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
            'database': 'reaxys_v2',
            'collection': 'buyables',
        }
        SOLVENTS = {
            'database': 'reaxys',
            'collection': 'solvents',
        }

        PREDICTOR = {
            'trained_model_path': '/home/ccoley/Make-It/makeit/predict/output/01_23_2017',
            'info': '01-23-17, model trained on 80k Reaxys examples, validated on 10k, tested on 10k. Nh1_200, Nh2_200, Nh3_200, l2_0, Nc_5000, enh_weight_0d1, context_weight_50, opt_adadelta, batch_5, moreFeatures'
        }

        database = db_client[RETRO_TRANSFORMS['database']]
        RETRO_DB = database[RETRO_TRANSFORMS['collection']]
        import makeit.retro.transformer as transformer

        RetroTransformer = transformer.Transformer(
            parallel=RETRO_TRANSFORMS['parallel'],
            nb_workers=RETRO_TRANSFORMS['nb_workers'],
        )
        print('Initialized retrotransformer')
        mincount_retro = RETRO_TRANSFORMS['mincount']
        RetroTransformer.load(RETRO_DB, mincount=mincount_retro, get_retro=True, get_synth=False)
        print('Loaded {} retro templates'.format(RetroTransformer.num_templates))
        RETRO_FOOTNOTE = 'Using {} retrosynthesis templates (mincount {}) from {}/{}'.format(RetroTransformer.num_templates,
                                                                                             mincount_retro,
                                                                                             RETRO_TRANSFORMS['database'],
                                                                                             RETRO_TRANSFORMS['collection'])

        ### Forward transformer
        database = db_client[SYNTH_TRANSFORMS['database']]
        SYNTH_DB = database[SYNTH_TRANSFORMS['collection']]
        SynthTransformer = transformer.Transformer()
        mincount_synth = SYNTH_TRANSFORMS['mincount']
        SynthTransformer.load(SYNTH_DB, mincount=100000000000000, get_retro=False, get_synth=True)
        print('Loaded {} forward templates'.format(SynthTransformer.num_templates))
        SYNTH_FOOTNOTE = 'Using {} forward templates (mincount {}) from {}/{}'.format(SynthTransformer.num_templates,
                                                                                      mincount_synth,
                                                                                      SYNTH_TRANSFORMS['database'],
                                                                                      SYNTH_TRANSFORMS['collection'])

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

        builder = TreeBuilder(Pricer=Pricer, RetroTransformer=RetroTransformer, nb_workers=4)

        # Intelligent predictor
        from makeit.webapp.forwardPredictor import ForwardPredictor

        predictor = ForwardPredictor(nb_workers=4, TRANSFORM_DB=SYNTH_DB, SOLVENT_DB=SOLVENT_DB)
        predictor.load_templates(mincount=mincount_synth)
        predictor.load_model(PREDICTOR['trained_model_path'])
        PREDICTOR_FOOTNOTE = 'Results generated using {} forward synthetic templates (mincount {}) from {}/{}, scored by a trained machine learning model: '.format(
            predictor.num_templates, mincount_synth, SYNTH_TRANSFORMS['database'], SYNTH_TRANSFORMS['collection']) + \
                             PREDICTOR['info']

        ## Context recommendation
        from askcos_site.functions.nnPredictor import NNConditionPredictor, lshf_nn, rxd_ids

        NN = NNConditionPredictor()
        NN.outputString = False

        import rdkit.Chem as Chem

    
    context_dict = {}
    plausibility_dict = {}
    plausible_trees = []
    done_trees = set()

    print('There are {} plausible tree(s) found'.format(len(plausible_trees)))
    attempt_number = 1
    while len(plausible_trees) < 3 and attempt_number <= int(args.max_expansions):

        FROOTROOT = os.path.join(os.path.dirname(__file__), 'CC')
        if not os.path.isdir(FROOTROOT):
            os.mkdir(FROOTROOT)
        FROOT = os.path.join(FROOTROOT, TARGET_LABEL)
        if not os.path.isdir(FROOT):
            os.mkdir(FROOT)
        fid = open(os.path.join(FROOT, '{}_integrated_log.txt'.format(attempt_number)), 'w')
        fid_success = open(os.path.join(FROOT, '{}_integrated_log_success.txt'.format(attempt_number)), 'w')
        fid_summary = open(os.path.join(FROOT, '{}_integrated_settings.txt'.format(attempt_number)), 'w')
        fid_summary.write('Target compound\t{}\n'.format(TARGET_LABEL))
        fid_summary.write('Target SMILES\t{}\n\n'.format(TARGET))
        fid_summary.write('Minimum retro template count\t{}\n'.format(args.retro_mincount))
        fid_summary.write('Minimum synth template count\t{}\n'.format(args.synth_mincount))
        fid_summary.write('Expansion time (s)\t{}\n'.format(expansion_time))
        fid_summary.write('Maximum depth\t{}\n'.format(max_depth))
        fid_summary.write('Maximum branching\t{}\n'.format(max_branching))
        fid_summary.write('Rank threshold for considering rxn plausible\t{}\n'.format(RANK_THRESHOLD_FOR_INCLUSION))
        fid_summary.write('Maximum trees considered\t{}\n\n'.format(max_trees))
        zero_time = time.time()

        def print_and_log(text, success=False):
            print(text)
            fid.write('[{%6.1f}] \t %s \n' % (time.time() - zero_time, text))
            if success:
                fid_success.write('[{%6.1f}] \t %s \n' % (time.time() - zero_time, text))
        
        print('This is expansion {}'.format(attempt_number))
        print('Current parameters: expansion_time={}, max_depth={}, max_branching={}, max_trees={}'.format(
            expansion_time, max_depth, max_branching, max_trees))

        main(TARGET=TARGET, expansion_time=expansion_time, max_depth=max_depth, max_branching=max_branching, max_trees=max_trees)

        fid_summary.write('Total processing time for this attempt (s)\t{}\n'.format(time.time() - zero_time))
        fid_summary.write('Total number of reactions evaluated so far\t{}\n'.format(len(plausibility_dict)))
        fid_summary.write('Total number of plausible trees found so far\t{}\n\n'.format(len(plausible_trees)))
        for plausible_tree in plausible_trees:
            fid_summary.write('{}\n\n'.format(plausible_tree))

        # Now update
        expansion_time *= 1.5 # Increasing time
        max_depth = min(max_depth + 2, 8) # 2 more steps up to 8
        max_branching += 5 # 5 more branches
        max_trees += 500 # 500 more trees to check
        attempt_number += 1

        # Save assessed-chemicals
        with open(os.path.join(FROOT, 'integrated_plog.txt'), 'w') as plog:
            plog.write('Target: {}.\nReaction smiles \t Context \t Plausibility \n'.format(TARGET))
            for rxn_smiles in plausibility_dict.keys():
                plog.write('{} \t {} \t {:.3f} \n'.format(
                    rxn_smiles, context_dict[rxn_smiles], plausibility_dict[rxn_smiles]
                ))

        fid.close()
        fid_success.close()
        fid_summary.close()

    quit(1)
