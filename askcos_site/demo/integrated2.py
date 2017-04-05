if __name__ == '__main__':

    from rdkit import RDLogger
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    import time
    import os
    fid = open(os.path.join(os.path.dirname(__file__), 'integrated_log.txt'), 'w')
    zero_time = time.time()

    def print_and_log(text):
        print(text)
        fid.write('[{%6.1f}] \t %s \n' % (time.time() - zero_time, text))

    # Atropine
    TARGET = 'CN3[C@H]1CC[C@@H]3C[C@@H](C1)OC(=O)C(CO)c2ccccc2'
    # # Fluconazole
    # TARGET = 'OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F'
    # Diphenhydramine
    #TARGET = 'O(CCN(C)C)C(c1ccccc1)c2ccccc2'


    expansion_time = 60
    max_depth = 5
    max_branching = 30
    max_trees = 1000


    from pymongo import MongoClient
    db_client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)

    RETRO_TRANSFORMS = {
        'database': 'reaxys',
        'collection': 'transforms_retro_v4', # 'lowe' or 'chematica'
        'mincount': 100,
        'parallel': False,
        'nb_workers': 1,
    }
    SYNTH_TRANSFORMS = {
        'database': 'reaxys',
        'collection': 'transforms_forward_v1'   ,
        'mincount': 50, 
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
    print_and_log('Initialized retrotransformer')
    mincount_retro = RETRO_TRANSFORMS['mincount']
    RetroTransformer.load(RETRO_DB, mincount = mincount_retro, get_retro = True, get_synth = False)
    print_and_log('Loaded {} retro templates'.format(RetroTransformer.num_templates))
    RETRO_FOOTNOTE = 'Using {} retrosynthesis templates (mincount {}) from {}/{}'.format(RetroTransformer.num_templates,
        mincount_retro, RETRO_TRANSFORMS['database'], RETRO_TRANSFORMS['collection'])

    ### Forward transformer 
    database = db_client[SYNTH_TRANSFORMS['database']]
    SYNTH_DB = database[SYNTH_TRANSFORMS['collection']]
    SynthTransformer = transformer.Transformer()
    mincount_synth = SYNTH_TRANSFORMS['mincount']
    SynthTransformer.load(SYNTH_DB, mincount = 100000000000000, get_retro = False, get_synth = True)
    print_and_log('Loaded {} forward templates'.format(SynthTransformer.num_templates))
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
    print_and_log('Loading prices...')
    import makeit.retro.pricer as pricer
    Pricer = pricer.Pricer()
    Pricer.load(CHEMICAL_DB, BUYABLE_DB)
    print_and_log('Loaded known prices')

    # Builder
    from makeit.webapp.treeBuilder import TreeBuilder 
    builder = TreeBuilder(Pricer = Pricer, RetroTransformer = RetroTransformer, nb_workers = 2)

    # Intelligent predictor
    from makeit.webapp.forwardPredictor import ForwardPredictor 
    predictor = ForwardPredictor(nb_workers = 2, TRANSFORM_DB = SYNTH_DB, SOLVENT_DB = SOLVENT_DB)
    predictor.load_templates(mincount = mincount_synth)
    predictor.load_model(PREDICTOR['trained_model_path'])
    PREDICTOR_FOOTNOTE = 'Results generated using {} forward synthetic templates (mincount {}) from {}/{}, scored by a trained machine learning model: '.format(predictor.num_templates,    mincount_synth, SYNTH_TRANSFORMS['database'], SYNTH_TRANSFORMS['collection']) + PREDICTOR['info']

    ## Context recommendation
    from askcos_site.functions.nnPredictor import nn_condition_predictor, lshf_nn, rxd_ids
    NN = nn_condition_predictor()
    NN.outputString = False

    import rdkit.Chem as Chem

    TARGET = Chem.MolToSmiles(Chem.MolFromSmiles(TARGET), isomericSmiles = False)
    print_and_log('Target: {}'.format(TARGET))
    builder.start_building(TARGET, max_depth = max_depth, max_branching = max_branching)
    print_and_log('Will wait {} seconds'.format(expansion_time))
    time.sleep(expansion_time)
    print_and_log(builder.info_string())

    builder.stop_building(timeout = 5)
    print_and_log('###########################################################')
    time.sleep(2)

    # Now define the reaction-checking function
    # Get context / evaluate?
    context_dict = {} 
    plausibility_dict = {}
    def check_this_reaction(tree, this_is_the_target=False):
        print_and_log('---------------')
        '''Evaluates this reaction and all its children'''
        global context_dict 
        global plausibility_dict 

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



            error = predictor.set_context(T = T1, reagents = rgt1, solvent = slvt1)
            if error is not None:
                print_and_log('Recommended context: {}'.format([T1, t1, y1, slvt1, rgt1, cat1]))
                print_and_log('NN recommended unrecognizable context for {}'.format(rxn_smiles))
                plausible = False

            else:
                # Run predictor
                predictor.run_in_foreground(reactants = '.'.join(reactants), intended_product = products[0], quit_if_unplausible = True)
                plausible = predictor.is_plausible()
                outcomes = predictor.return_top()

                print_and_log('necessary reagents: {}'.format(necessary_reagent))
                print_and_log('Recommended context: {}'.format([T1, t1, y1, slvt1, rgt1, cat1]))
                print_and_log('plausible? {}'.format(plausible))
                print_and_log(outcomes[:3])

            if plausible: 
                for outcome in outcomes:
                    if outcome['smiles'] == products[0]:
                        plausible = float(outcome['prob'])
                        break
                plausibility_dict[rxn_smiles] = plausible # save probability directly

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
                print_and_log('################################################')
                print_and_log('Found completely plausible tree!')
                print_and_log(tree)
                print_and_log('################################################')
            elif this_is_the_target:
                print_and_log('Found tree with unfeasible children...thats not helpful')
        elif this_is_the_target:
            print_and_log('Found tree with unfeasible first step...thats not helpful')


        return plausible

    # Define chemical and reaction helper functions (redundant with builder)
    def chem_dict(_id, smiles, ppg, children = []):
        '''Chemical object as expected by website'''
        return {
            'id': _id,
            'is_chemical': True,
            'smiles' : smiles,
            'ppg' : ppg,
            'children': children,
        }

    def rxn_dict(_id, info, necessary_reagent = '', num_examples = 0, children = []):
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
                    yield [] # Nope, this is a viable node
                else:
                    # Try going deeper via DLS_rxn function
                    for rxn_id in tree_dict[chem_id]['prod_of']:
                        rxn_info_string = ''
                        for path in DLS_rxn(rxn_id, depth):
                            yield [rxn_dict(rxn_id, rxn_info_string, tree_dict[rxn_id]['necessary_reagent'], tree_dict[rxn_id]['num_examples'], children = path)]

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
                        chem_dict(chem_id, tree_dict[chem_id]['smiles'], tree_dict[chem_id]['ppg'], children = path)
                    ]

            # Two reactants? want to capture all combinations of each node's options
            elif len(generators) == 2:
                chem_id0 = tree_dict[rxn_id]['rcts'][0]
                chem_id1 = tree_dict[rxn_id]['rcts'][1]
                for path0 in generators[0]:
                    for path1 in generators[1]:
                        yield [
                            chem_dict(chem_id0, tree_dict[chem_id0]['smiles'], tree_dict[chem_id0]['ppg'], children = path0),
                            chem_dict(chem_id1, tree_dict[chem_id1]['smiles'], tree_dict[chem_id1]['ppg'], children = path1)
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
                                chem_dict(chem_id0, tree_dict[chem_id0]['smiles'], tree_dict[chem_id0]['ppg'], children = path0),
                                chem_dict(chem_id1, tree_dict[chem_id1]['smiles'], tree_dict[chem_id1]['ppg'], children = path1),
                                chem_dict(chem_id2, tree_dict[chem_id2]['smiles'], tree_dict[chem_id2]['ppg'], children = path2),
                            ]

            else:
                print_and_log('Too many reactants! Only have cases 1-3 programmed')
                raise ValueError('Too many reactants! Only have cases 1-3 programmed')

        for depth in range(max_depth):
            print_and_log('~~ looking at depth {}'.format(depth))
            for path in DLS_chem(1, depth, headNode=True):
                yield chem_dict(1, tree_dict[1]['smiles'], tree_dict[1]['ppg'], children=path)

    # Look for trees
    print_and_log('Looking for trees now')
    import hashlib 
    import json
    done_trees = set()
    trees = []
    counter = 0
    for tree in IDDFS_generator(
            tree_dict=builder.tree_dict, 
            max_depth=builder.max_depth):
        # Step 1 - check if unique
        hashkey = hashlib.sha1(json.dumps(tree, sort_keys=True)).hexdigest()
        if hashkey in done_trees:
            print_and_log('Found duplicate tree, skipping')
            continue
        # Step 2 - check feasibility
        check_this_reaction(tree, this_is_the_target=True)

        counter = counter + 1
        if counter == max_trees: 
            break

    print_and_log('Got out of final IDDFS_generator loop, done with program')
    builder.stop_building()
    predictor.stop()
    print_and_log('Stopped builder and predictor')
    del builder 
    del predictor
    del NN
    fid.close()
    quit(1)