import json
from rdkit import Chem
from collections import defaultdict
from django.http import JsonResponse
from askcos_site.askcos_celery.treebuilder.tb_coordinator_mcts import get_buyable_paths as get_buyable_paths_mcts

def tree_builder(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    run_async = request.GET.get('async', True)
    orig_smiles = request.GET.get('smiles', None)
    mol = Chem.MolFromSmiles(orig_smiles)
    smiles = Chem.MolToSmiles(mol)
    max_depth = int(request.GET.get('max_depth', 4))
    max_branching = int(request.GET.get('max_branching', 25))
    expansion_time = int(request.GET.get('expansion_time', 60))
    max_ppg = int(request.GET.get('max_ppg', 10))
    template_count = int(request.GET.get('template_count', '100'))
    max_cum_prob = float(request.GET.get('max_cum_prob', '0.995'))
    chemical_property_logic = str(request.GET.get('chemical_property_logic', 'none'))
    max_chemprop_c = int(request.GET.get('max_chemprop_c', '0'))
    max_chemprop_n = int(request.GET.get('max_chemprop_n', '0'))
    max_chemprop_o = int(request.GET.get('max_chemprop_o', '0'))
    max_chemprop_h = int(request.GET.get('max_chemprop_h', '0'))
    chemical_popularity_logic = str(request.GET.get('chemical_popularity_logic', 'none'))
    min_chempop_reactants = int(request.GET.get('min_chempop_reactants', 5))
    min_chempop_products = int(request.GET.get('min_chempop_products', 5))
    filter_threshold = float(request.GET.get('filter_threshold', 0.75))
    apply_fast_filter = filter_threshold > 0
    return_first = json.loads(request.GET.get('return_first', 'false'))
    
    blacklisted_reactions = request.GET.get('blacklisted_reactions', [])
    forbidden_molecules = request.GET.get('forbidden_molecules', [])
    
    default_val = 1e9 if chemical_property_logic == 'and' else 0
    max_natom_dict = defaultdict(lambda: default_val, {
        'logic': chemical_property_logic,
        'C': max_chemprop_c,
        'N': max_chemprop_n,
        'O': max_chemprop_o,
        'H': max_chemprop_h,
    })
    min_chemical_history_dict = {
        'logic': chemical_popularity_logic,
        'as_reactant': min_chempop_reactants,
        'as_product': min_chempop_products,
    }
    
    res = get_buyable_paths_mcts.delay(smiles, max_branching=max_branching, max_depth=max_depth,
                                  max_ppg=max_ppg, expansion_time=expansion_time, max_trees=500,
                                  known_bad_reactions=blacklisted_reactions,
                                  forbidden_molecules=forbidden_molecules,
                                  max_cum_template_prob=max_cum_prob, template_count=template_count,
                                  max_natom_dict=max_natom_dict, min_chemical_history_dict=min_chemical_history_dict,
                                  apply_fast_filter=apply_fast_filter, filter_threshold=filter_threshold,
                                  return_first=return_first)
    
    if run_async:
        resp['id'] = res.id
        resp['state'] = res.state
        return JsonResponse(resp)
    
    try:
        (tree_status, trees) = res.get(expansion_time * 3)
    except:
        resp['error'] = 'API request timed out (after {})'.format(expansion_time * 3)
        res.revoke()
        return JsonResponse(resp)
    
    resp['trees'] = trees
    
    return JsonResponse(resp)