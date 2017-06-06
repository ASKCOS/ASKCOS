import os
import sys
import time
import argparse
import numpy as np
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
import rdkit.Chem as Chem 

def wait_to_get_result(res, timeout=120, update_freq=5, sleep=0.1, 
    label='result', default=None):
    '''Helper func to wrap around a .get call'''
    time_start = time.time()
    time_last_update = time_start - update_freq
    try:
        while True:
            time_now = time.time()
            if res.ready():
                return res.get()
            if time_now - time_start > timeout:
                res.revoke(terminate=True)
                print('Result for {} took too long!'.format(label))
                return default
            if time_now - time_last_update >= update_freq:
                print('Waiting for {}, {:.2f}s elapsed, timeout at {}'.format(
                    label, time_now - time_start, timeout
                ))
                time_last_update = time_now
            time.sleep(sleep)
    except KeyboardInterrupt:
        res.revoke(terminate=True)
        # res.forget() # gracefully(?) exit
        raise(KeyboardInterrupt('Propagating interrupt after revoking most recent task!'))


def main(TARGET, expansion_time=60.0, max_depth=5, max_branching=30, max_trees=1000,
        max_ppg=1e8):
    global context_dict
    global plausibility_dict
    global plausible_trees
    global print_and_log
    global log_plausibility

    # Function to recontruct the plausible tree
    def construct_plausible_tree(tree, depth=0):
        overall_plausibility = 1.0
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
            overall_plausibility *= plausibility_dict[rxn_smiles]
            depth += 1
            for child in rxn['children']:
                overall_plausibility *= construct_plausible_tree(child, depth=depth)
        return overall_plausibility

    # Helper
    def fix_rgt_cat_slvt(rgt1, cat1, slvt1):
        # Merge cat and reagent for forward predictor
        if rgt1 and cat1:
            rgt1 = rgt1 + '.' + cat1
        elif cat1:
            rgt1 = cat1
        # Reduce solvent to single one
        if '.' in slvt1:
            slvt1 = slvt1.split('.')[0]
        return (rgt1, cat1, slvt1)

    def trim_trailing_period(txt):
        if txt:
            if txt[-1] == '.':
                return txt[:-1]
        return txt

    def check_this_reaction(reactants, products, necessary_reagent):
        rxn_smiles = '.'.join(sorted(reactants)) + '>>' + products[0]
        try:
            return plausibility_dict[rxn_smiles]
        except KeyError:
            pass 

        # It is easier to test multiple contexts when there is no necessary_reagent
        if not necessary_reagent: 
            print_and_log('no necessary reagent needed - using multicontext approach')
            contexts_res = get_context_recommendation.delay(rxn=[reactants, products], n=MAX_NUM_CONTEXTS)
            contexts = contexts_res.get(30)
            contexts_for_predictor = []
            for (T1, t1, y1, slvt1, rgt1, cat1) in contexts:
                slvt1 = trim_trailing_period(slvt1)
                rgt1 = trim_trailing_period(rgt1)
                cat1 = trim_trailing_period(cat1)
                (rgt1, cat1, slvt1) = fix_rgt_cat_slvt(rgt1, cat1, slvt1)
                contexts_for_predictor.append((T1, rgt1, slvt1))

            print('Starting forward predictor')
            all_outcomes_res = get_outcomes.delay('.'.join(reactants), contexts_for_predictor,
                top_n=RANK_THRESHOLD_FOR_INCLUSION)
            all_outcomes = wait_to_get_result(all_outcomes_res, timeout=300, update_freq=5, sleep=0.1, 
                label='forward predictor results (multi-context)',
                default=[[] for i in range(len(contexts_for_predictor))])
            print('Done!')
            all_outcomes = all_outcomes_res.get()
            plausible = [0. for i in range(len(all_outcomes))]
            for i, outcomes in enumerate(all_outcomes):
                for outcome in outcomes:
                    if outcome['smiles'] == products[0]:
                        plausible[i] = float(outcome['prob'])
                        break
            best_context_i = np.argmax(plausible)
            plausible = plausible[best_context_i]
            context_dict[rxn_smiles] = contexts_for_predictor[best_context_i]
            print_and_log('Recommended BEST context out of {} attempted: {}'.format(len(contexts), contexts[best_context_i]))
            print_and_log('plausible? {}'.format(plausible))
            print_and_log(all_outcomes[best_context_i][:min(3, len(all_outcomes[best_context_i]))])

        else:
            print_and_log('necessary reagent: {}'.format(necessary_reagent))
            context_res = get_context_recommendation.delay(rxn=[reactants, products], n=1)
            [T1, t1, y1, slvt1, rgt1, cat1] = context_res.get(30)
            slvt1 = trim_trailing_period(slvt1)
            rgt1 = trim_trailing_period(rgt1)
            cat1 = trim_trailing_period(cat1)
            context_dict[rxn_smiles] = [T1, t1, y1, slvt1, rgt1, cat1]

            # Add reagent to reactants for the purpose of forward prediction
            # this is important for necessary reagents, e.g., chlorine source
            if not rgt1 and necessary_reagent:
                #TODO: add a match to see if the rgt1 is even capable of providing the right atoms
                print_and_log('necessary reagent not found!')
                plausibility_dict[rxn_smiles] = 0.
                return 0.
            if rgt1:
                reactants.append(rgt1)
                print('Added reagent {} to reactant pool'.format(rgt1))
            
            (rgt1, cat1, slvt1) = fix_rgt_cat_slvt(rgt1, cat1, slvt1)
            
            print('Starting forward predictor (one context)')
            all_outcomes_res = get_outcomes.delay('.'.join(reactants), 
                [(T1, rgt1, slvt1)], top_n=RANK_THRESHOLD_FOR_INCLUSION)
            all_outcomes = wait_to_get_result(all_outcomes_res, timeout=120, update_freq=5, sleep=0.1, 
                label='forward predictor results (single-context)',
                default=[[]])
            print('Done!')
            outcomes = all_outcomes[0]
            plausible = 0.
            for outcome in outcomes:
                if outcome['smiles'] == products[0]:
                    plausible = float(outcome['prob'])
                    break
            print_and_log('necessary reagents: {}'.format(necessary_reagent))
            print_and_log('Recommended context: {}'.format([T1, t1, y1, slvt1, rgt1, cat1]))
            print_and_log('plausible? {}'.format(plausible))
            print_and_log(outcomes[:min(3, len(outcomes))])

        if plausible < PROB_THRESHOLD_FOR_INCLUSION:
            plausible = 0.
        plausibility_dict[rxn_smiles] = plausible 
        log_plausibility(rxn_smiles)
        
        return plausible

    # Now define the reaction-checking function
    # Get context / evaluate?
    def check_this_tree(tree, this_is_the_target=False):
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

        plausible = check_this_reaction(reactants, products, necessary_reagent)

        # Look at children?
        if plausible > 0:
            all_plausible = plausible
            for child in rxn['children']:
                all_plausible *= check_this_tree(child)

            # Now report
            if all_plausible and this_is_the_target:
                plausible_trees.append(tree)
                print_and_log('################################################', success=True)
                print_and_log('Found completely plausible tree!', success=True)
                print_and_log('Overall plausibility for whole tree: {}'.format(all_plausible), success=True)
                print_and_log(tree, success=True)
                construct_plausible_tree(tree, depth=0)
                print_and_log('################################################', success=True)
            elif this_is_the_target:
                print_and_log('Found tree with unfeasible children...thats not helpful')

            return all_plausible # need to propagate up

        elif this_is_the_target:
            print_and_log('Found tree with unfeasible first step...thats not helpful')

        return plausible

    print('Starting expansion of tree')
    known_bad_reactions = [key for (key, val) in plausibility_dict.iteritems() if val == 0.]
    tree_result_async = get_buyable_paths.delay(TARGET, max_time=expansion_time, max_depth=max_depth, 
        max_branching=max_branching, max_trees=max_trees, known_bad_reactions=known_bad_reactions,
        return_d1_if_no_trees=True, max_ppg=max_ppg)
    tree_result = wait_to_get_result(tree_result_async, timeout=expansion_time*10., update_freq=5,
        default=((0,0,{}),[],[]),
        label='trees for target {} (t={:.1f}, d={}, b={}, max_trees={}, known_bad_reactions={})'.format(
            TARGET, expansion_time, max_depth, max_branching, max_trees, len(known_bad_reactions)
        ))

    # NORMAL SITUATION: some pathways found!
    d1_reactions = None
    if len(tree_result) == 2:
        (tree_status, trees) = tree_result 
    elif len(tree_result) == 3:
        (tree_status, trees, d1_reactions) = tree_result 
    else:
        raise ValueError('result from tree builder has unexpected size')

    # Save status
    (num_chemicals, num_reactions, at_depth) = tree_status
    print_and_log('After expanding, {} total chemicals and {} total reactions'.format(num_chemicals, num_reactions),
        success=True, summary=True)
    for (depth, count) in sorted(at_depth.iteritems(), key=lambda x: x[0]):
        print_and_log('At depth {}, {} nodes'.format(depth, count), success=True, summary=True)

    # Evaluate
    for tree in trees:
        # Step 1 - check if with exception
        if tree is None:
            print_and_log('Found tree with exceptions, skipping')
            continue
        # Step 2 - check feasibility
        check_this_tree(tree, this_is_the_target=True)
    print_and_log('Done with all trees found')

    # Evaluate NON-trees
    if d1_reactions is not None:
        print_and_log('Because no trees found, will look at d=1 reactions', success=True)
        for (reactants, products, necessary_reagent) in d1_reactions:
            rxn_smiles = '.'.join(sorted(reactants)) + '>>' + products[0]
            print_and_log('reaction smiles: {}'.format(rxn_smiles))
            plausible = check_this_reaction(reactants, products, necessary_reagent)
            if plausible:
                print_and_log('Found one plausible d=1 reaction - good enough to continue')
                break # one good one is enough to try...

    return

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
    parser.add_argument('--max_trees', type=int, default=10000,
                        help='Maxinum number of trees to examine, default 10000')
    parser.add_argument('--retro_mincount', type=int, default=100,
                        help = 'Minimum template count for retrosynthetic templates, default 100')
    parser.add_argument('--synth_mincount', type=int, default=50,
                        help = 'Minimum template count for synthetic templates, default 50')
    parser.add_argument('--rank_threshold', type=int, default=10,
                        help = 'Rank threshold for considering a reaction to be plausible, default 10')
    parser.add_argument('--p_threshold', type=float, default=0.025,
                        help = 'Probability threshold for considering a reaction to be plausible, default 0.025')
    parser.add_argument('--max_expansions', type=int, default=8,
                        help = 'Maximum number of search expansions to try, default 8')
    parser.add_argument('--max_contexts', type=int, default=10,
                        help = 'Maximum number of contexts to try when reagents not necessary, default 10')
    parser.add_argument('--max_ppg', type=int, default=10,
                        help = 'Maximum price per gram for a buyable chemical, default 10')
    parser.add_argument('--min_trees_success', type=int, default=5,
                        help = 'Minimum number of trees before quitting, default 5')
    args = parser.parse_args()

    

    expansion_time = float(args.expansion_time)
    max_depth = int(args.max_depth)
    max_branching = int(args.max_branching)
    max_trees = int(args.max_trees)
    TARGET = str(args.TARGET)
    RANK_THRESHOLD_FOR_INCLUSION = int(args.rank_threshold)
    PROB_THRESHOLD_FOR_INCLUSION = float(args.p_threshold)
    MAX_NUM_CONTEXTS = int(args.max_contexts)
    max_ppg = int(args.max_ppg)
    min_trees_success = int(args.min_trees_success)

    import urllib2
    smiles = urllib2.urlopen('https://cactus.nci.nih.gov/chemical/structure/{}/smiles'.format(TARGET)).read()
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError('Could not resolve SMILES ({}) in a way parseable by RDKit'.format(smiles))
    smiles = Chem.MolToSmiles(mol)
    print('Resolved {} smiles -> {}'.format(TARGET, smiles))
    TARGET_LABEL = TARGET 
    TARGET = smiles

    FROOTROOT = os.path.join(os.path.dirname(__file__), 'CC')
    if not os.path.isdir(FROOTROOT):
        os.mkdir(FROOTROOT)
    FROOT = os.path.join(FROOTROOT, TARGET_LABEL)
    if not os.path.isdir(FROOT):
        os.mkdir(FROOT)

    from askcos_site.askcos_celery.contextrecommender.worker import get_context_recommendation
    res = get_context_recommendation.delay([['CCO','CCBr'],['CCOCC']], n=1)
    wait_to_get_result(res, timeout=120, update_freq=1, sleep=0.1, 
        label='context recommender to come online')

    from askcos_site.askcos_celery.treebuilder.worker import get_top_precursors
    res = get_top_precursors.delay('CCCCCOCCCNC(=O)CCC')
    wait_to_get_result(res, timeout=120, update_freq=1, sleep=0.1, 
        label='a tree builder worker to come online')

    from askcos_site.askcos_celery.treebuilder.coordinator import get_buyable_paths
    res = get_buyable_paths.delay('CCCOCCC', mincount=0, max_branching=10, max_depth=2, 
        max_ppg=1e8, max_time=10, max_trees=3, reporting_freq=5,
        known_bad_reactions=['CCCO>>CCCOCCC'])
    wait_to_get_result(res, timeout=120, update_freq=1, sleep=0.1, 
        label='a tree builder coordinator to come online')

    from askcos_site.askcos_celery.forwardpredictor.worker import get_candidate_edits
    res = get_candidate_edits.delay(reactants_smiles='[C:1][C:2][O:3].[Br:4][C:5]', start_at=0, end_at=100)
    wait_to_get_result(res, timeout=120, update_freq=1, sleep=0.1, 
        label='a forward predictor worker to come online')

    from askcos_site.askcos_celery.forwardpredictor.coordinator import get_outcomes
    res = get_outcomes.delay('CCCO.CCCBr', contexts=[('20','[Na+].[OH-]', 'water'), ('30', '', 'toluene')], mincount=100)
    wait_to_get_result(res, timeout=120, update_freq=1, sleep=0.1, 
        label='a forward predictor coordinator to come online')

    context_dict = {}
    plausibility_dict = {}
    plausible_trees = []
    done_trees = set()

    plog = open(os.path.join(FROOT, 'integrated_plog.txt'), 'w')
    plog.write('Target: {} ({})'.format(TARGET_LABEL, TARGET))
    plog.write('{}\t{}\t{}\n'.format(
        'Reaction SMILES', 'Context', 'Plausibility'
    ))
    def log_plausibility(rxn_smiles):
        try:
            plog.write('{}\t{}\t{:.4f}\n'.format(
                rxn_smiles, context_dict[rxn_smiles], plausibility_dict[rxn_smiles]
            ))
        except KeyError:
            print('Cannot log reaction smiles {}'.format(rxn_smiles))

    print('There are {} plausible tree(s) found'.format(len(plausible_trees)))
    attempt_number = 1
    while len(plausible_trees) < min_trees_success and attempt_number <= int(args.max_expansions):

        fid = open(os.path.join(FROOT, '{}_integrated_log.txt'.format(attempt_number)), 'w')
        fid_success = open(os.path.join(FROOT, '{}_integrated_log_success.txt'.format(attempt_number)), 'w')
        fid_summary = open(os.path.join(FROOT, '{}_integrated_summary.txt'.format(attempt_number)), 'w')
        fid_summary.write('Target compound\t{}\n'.format(TARGET_LABEL))
        fid_summary.write('Target SMILES\t{}\n\n'.format(TARGET))
        fid_summary.write('Minimum retro template count\t{}\n'.format(args.retro_mincount))
        fid_summary.write('Minimum synth template count\t{}\n'.format(args.synth_mincount))
        fid_summary.write('Expansion time (s)\t{}\n'.format(expansion_time))
        fid_summary.write('Maximum depth\t{}\n'.format(max_depth))
        fid_summary.write('Maximum branching\t{}\n'.format(max_branching))
        fid_summary.write('Rank threshold for considering rxn plausible\t{}\n'.format(RANK_THRESHOLD_FOR_INCLUSION))
        fid_summary.write('Probability threshold for considering rxn plausible\t{}\n'.format(PROB_THRESHOLD_FOR_INCLUSION))
        fid_summary.write('Number of contexts recommended when no necessary_reagent\t{}\n'.format(MAX_NUM_CONTEXTS))
        fid_summary.write('Maximum trees considered\t{}\n'.format(max_trees))
        fid_summary.write('Maximum ppg for buyability\t{}\n'.format(max_ppg))
        fid_summary.write('\n')
        zero_time = time.time()

        def print_and_log(text, success=False, summary=False):
            print(text)
            fid.write('[{%6.1f}] \t %s \n' % (time.time() - zero_time, text))
            if success:
                fid_success.write('[{%6.1f}] \t %s \n' % (time.time() - zero_time, text))
            if summary:
                fid_summary.write('{}\n'.format(text))

        print('This is expansion {}'.format(attempt_number))
        print('Current parameters: expansion_time={}, max_depth={}, max_branching={}, max_trees={}'.format(
            expansion_time, max_depth, max_branching, max_trees))

        main(TARGET=TARGET, expansion_time=expansion_time, max_depth=max_depth, 
            max_branching=max_branching, max_trees=max_trees, max_ppg=max_ppg)

        fid_summary.write('\nTotal processing time for this attempt (s)\t{}\n'.format(time.time() - zero_time))
        fid_summary.write('Total number of reactions evaluated so far\t{}\n'.format(len(plausibility_dict)))
        fid_summary.write('Total number of plausible trees found so far\t{}\n\n'.format(len(plausible_trees)))
        for plausible_tree in plausible_trees:
            fid_summary.write('{}\n\n'.format(plausible_tree))

        # Now update
        if attempt_number % 2 == 1:
            expansion_time *= 1.5 # Increasing time
        else:
            max_depth = min(max_depth + 1, 8) # 2 more steps up to 8
        max_branching += 5 # 5 more branches
        max_trees += 5000 # 5000 more trees to check - shouldn't ever be limiting
        attempt_number += 1

        fid.close()
        fid_success.close()
        fid_summary.close()
    plog.close()

    quit(1)
