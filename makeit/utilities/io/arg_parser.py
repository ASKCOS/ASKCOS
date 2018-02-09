import argparse
import makeit.global_config as gc
arg_parser_loc = 'arg_parser'


def setup_parser():

    parser = argparse.ArgumentParser()

    # Initialize default values for the different arguments
    expansion_time = 60.0
    max_depth = 4
    max_branching = 20
    max_trees = 5000
    context_recommender = gc.neural_network
    forward_scoring_method = gc.templatefree
    tree_scoring_method = gc.forwardonly
    context_prioritization = gc.probability
    template_prioritization = gc.relevance
    precursor_prioritization = gc.heuristic
    retro_mincount = 25
    retro_mincount_chiral = 10
    synth_mincount = 25
    rank_threshold = 5
    p_threshold = 0.01
    max_expansions = 8
    max_contexts = 10
    max_ppg = 100
    min_trees_success = 5
    precursor_score_mode = gc.max
    max_cum_template_prob = 0.995 
    output_dir = 'output'
    chiral = True
    nproc = 2
    celery = False
    template_count = 100
    parallel_tree = False

    # Set all arguments
    parser.add_argument('--TARGET', type=str,
                        help='Name of the target chemical.')
    parser.add_argument('--expansion_time', type=float, default=expansion_time,
                        help='Retrosynthesis tree expansion time. Default expansion time is {} s.'.format(expansion_time))
    parser.add_argument('--max_depth', type=int, default=max_depth,
                        help='Maximum depth for retrosynthesis tree expansion. Default value is {}.'.format(max_depth))
    parser.add_argument('--max_branching', type=int, default=max_branching,
                        help='Maximal branching factor for retrosynthesis tree expansion. Default value is {}.'.format(max_branching))
    parser.add_argument('--max_trees', type=int, default=max_trees,
                        help='Maximum number of trees that will be investigated. Default value is {}.'.format(max_trees))
    parser.add_argument('--retro_mincount', type=int, default=retro_mincount,
                        help='Minimum template count for retrosynthetic templates. Default value is {}.'.format(retro_mincount))
    parser.add_argument('--retro_mincount_chiral', type=int, default=retro_mincount_chiral,
                        help='Minimum template count for chiral retrosynthetic templates. Default value is {}'.format(retro_mincount_chiral))
    parser.add_argument('--synth_mincount', type=int,  default=synth_mincount,
                        help='Minimum template count for synthetic templates. Default value is {}.'.format(synth_mincount))
    parser.add_argument('--rank_threshold', type=int, default=rank_threshold,
                        help='Maximal rank for which a reaction is still considered plausible. Default value is {}.'.format(rank_threshold))
    parser.add_argument('--prob_threshold', type=int, default=p_threshold,
                        help='Minimal probability for which a reaction is still considered to be plausible. Default value is {}.'.format(p_threshold))
    parser.add_argument('--max_expansions', type=int, default=max_expansions,
                        help='Maximum number of search expansions to attempt. Default value is {}.'.format(max_expansions))
    parser.add_argument('--max_contexts', type=int, default=max_contexts,
                        help='Maximum number of contexts to try when reagents are not necessary. Default value is {}.'.format(max_contexts))
    parser.add_argument('--template_count', type=int, default=template_count,
                        help='Maximum number of templates to be applied for each expansion step. Default value is {}'.format(template_count))
    parser.add_argument('--max_ppg', type=int, default=max_ppg,
                        help='Maximum price per gram to consider a chemical buyable. Default value is {}.'.format(max_ppg))
    parser.add_argument('--min_trees_success', type=int, default=min_trees_success,
                        help='Minimum number of trees to evaluate before quitting. Default value is {}.'.format(min_trees_success))
    parser.add_argument('--output', type=str, default=output_dir,
                        help='Specify the directory to which all output should be written. Default directory is {}'.format(output_dir))
    parser.add_argument('--chiral', type=bool, default=chiral,
                        help='Specify whether a chiral transformer should be used or not. Default is {}'.format(chiral))
    parser.add_argument('--nproc', type=int, default=nproc,
                        help='Specify the number of processors on which the case should be parallelized. Default is single-core (1).')
    parser.add_argument('--celery', type=bool, default=celery,
                        help='Specify whether preloaded celery worker nodes should be used. Default is {}'.format(celery))
    parser.add_argument('--context_recommender', type=str, default=context_recommender,
                        help='Specify which method for context recommendation should be used. Default is {}'.format(context_recommender))
    parser.add_argument('--context_prioritization', type=str, default=context_prioritization,
                        help='Specify by which parameter the top context should be determined. Default is {}'.format(context_prioritization))
    parser.add_argument('--forward_scoring', type=str, default=forward_scoring_method,
                        help='Specify which method is to be used to evaluate a forward reaction step. Default is {}'.format(forward_scoring_method))
    parser.add_argument('--tree_scoring', type=str, default=tree_scoring_method,
                        help='Specify how the score of a tree should be determined. Default is {}'.format(tree_scoring_method))
    parser.add_argument('--template_prioritization', type=str, default=template_prioritization,
                        help='Specify which method should be used to order templates. Default is {}'.format(template_prioritization))
    parser.add_argument('--precursor_prioritization', type=str, default=precursor_prioritization,
                        help='Specify which method should be used to order the found precursors. Default is {}'.format(precursor_prioritization))
    parser.add_argument('--parallel_tree', type=bool, default=parallel_tree,
                        help='Specify whether tree parallelization should be carried out in parallel. Default is{}'.format(parallel_tree))
    parser.add_argument('--precursor_score_mode', type=str, default=precursor_score_mode,
                        help='Specify how the score of a list of precursors should be determined when using synthetic complexity. Default is {}.'.format(precursor_score_mode))
    parser.add_argument('--max_cum_template_prob', type=float, default = max_cum_template_prob,
                        help='Specify the cumulative template probability cut-off. Default is {}'.format(max_cum_template_prob))
    return parser


def get_args():

    parser = setup_parser()
    args = parser.parse_args()

    return args
