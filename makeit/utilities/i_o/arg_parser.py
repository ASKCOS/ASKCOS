import argparse
arg_parser_loc = 'arg_parser'
def setup_parser():
    
    parser = argparse.ArgumentParser()
    
    #Initialize default values for the different arguments
    expansion_time = 60.0
    max_depth = 3
    max_branching = 10
    max_trees = 10000
    retro_mincount = 100
    synth_mincount = 50
    rank_threshold = 10
    p_threshold = 0.025
    max_expansions = 8
    max_contexts = 10
    max_ppg = 10
    min_trees_success = 5
    output_dir = 'output'
    chiral = False
    nproc = 1
    celery = False
    
    
    #Set all arguments
    parser.add_argument('--TARGET', type = str, 
                        help = 'Name of the target chemical.')
    parser.add_argument('--expansion_time', type = float, default = expansion_time, 
                        help = 'Retrosynthesis tree expansion time. Default expansion time is {} s.'.format(expansion_time))
    parser.add_argument('--max_depth', type = int, default = max_depth,
                        help = 'Maximum depth for retrosynthesis tree expansion. Default value is {}.'.format(max_depth))
    parser.add_argument('--max_branching', type = int, default = max_branching,
                        help = 'Maximal branching factor for retrosynthesis tree expansion. Default value is {}.'.format(max_branching))
    parser.add_argument('--max_trees', type = int, default = max_trees,
                        help = 'Maximum number of trees that will be investigated. Default value is {}.'.format(max_trees))
    parser.add_argument('--retro_mincount', type = int, default = retro_mincount,
                        help = 'Minimum template count for retrosynthetic templates. Default value is {}.'.format(retro_mincount))
    parser.add_argument('--synth_mincount', type = int,  default = synth_mincount,
                        help = 'Minimum template count for synthetic templates. Default value is {}.'.format(synth_mincount))
    parser.add_argument('--rank_threshold', type = int, default = rank_threshold,
                        help = 'Maximal rank for which a reaction is still considered plausible. Default value is {}.'.format(rank_threshold))
    parser.add_argument('--p_threshold', type = int, default = p_threshold,
                        help = 'Minimal probability for which a reaction is still considered to be plausible. Default value is {}.'.format(p_threshold))
    parser.add_argument('--max_expansions', type = int, default = max_expansions,
                        help = 'Maximum number of search expansions to attempt. Default value is {}.'.format(max_expansions))
    parser.add_argument('--max_contexts', type = int, default = max_contexts,
                        help = 'Maximum number of contexts to try when reagents are not necessary. Default value is {}.'.format(max_contexts))
    parser.add_argument('--max_ppg', type = int, default = max_ppg,
                        help = 'Maximum price per gram to consider a chemical buyable. Default value is {}.'.format(max_ppg))
    parser.add_argument('--min_trees_success', type = int, default = min_trees_success,
                        help = 'Minimum number of trees to evaluate before quitting. Default value is {}.'.format(min_trees_success))
    parser.add_argument('--output', type = str, default = output_dir,
                        help = 'Specify the directory to which all output should be written. Default directory is {}'.format(output_dir))
    parser.add_argument('--chiral', type = boolean, default = chiral,
                        help = 'Specify whether a chiral transformer should be used or not. Default is {}'.format(chiral))
    parser.add_argument('--nproc', type = int, default = nproc,
                        help = 'Specify the number of processors on which the case should be parallelized. Default is single-core (1).')
    parser.add_argument('--celery', type = boolean, default = celery,
                        help = 'Specify whether preloaded celery worker nodes should be used. Default is {}'.format(celery))
    return parser

def get_args():
    
    parser = setup_parser()
    args = parser.parse_args()
    
    return args