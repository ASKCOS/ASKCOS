

```python
import requests
from pprint import pprint
```


```python
HOST = 'http://<your-askcos-ip>'
```

# ASKCOS Web API

## Endpoints

API endpoints exist for the following services (described in more depth below):

|Endpoint|Description|
|:---|---:|
|/api/retro/|Single-step retrosynthetic prediction for a given target molecule|
|/api/context/|Context recommender for a reaction|
|/api/forward/|Reaction (product) prediction given reactants and context|
|/api/treebuilder/|Synthetic pathway tree builder|
|/api/template/|Retrosynthetic template lookup|
|/api/fast-filter/|Coarse-filter for reactions plausibility|
|/api/scscore/|Synthetic complexity score|
|/api/price/|Buyable price|
|/api/celery/|Query the status of the celery workers|

## Making requests

GET requests can be made passing in required and optional parameters (service-specific parameters are described below for each endpoint). In many cases, you will be providing SMILES strings, which often may not play well with HTTP. You should always escape your parameters. In the examples below, the `requests` python package allows you to provide the parameters as a dictionary of key-value pairs and takes care of the url escaping for you. If you create your own url queries, please remember to escape your SMILES strings.

In the examples below, a dictionary named `params` will be built to contain the parameters for each service. `requests.get` will be used, giving the API endpoint and the params dictionary. This returns a `Response` object, and the final results can be obtained using `Response.json()`. The `pprint` module is used to nicely format the JSON results.

For each `params` dictionary, required key-value pairs will be identified with `# required`, and the rest of the values can be considered optional, as the values used here are the default values used (unless otherwise stated). In these examples, the number of results being returned is being significantly decreased to make this document more readable.

## /api/retro/  
Given a `target` SMILES string, predict precursors. The optional settings adjust how precursors are generated, or filtered after generation.  
  
The response will be a json with the original request variables as well as `precursors`, a list of the suggested results. The actualsuggested precursor smiles strings are stored in `smiles` (molecules concatenated with '.'), and `smiles_split` (list of molecules).


```python
params = {
    # required
    'target': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
    
    # optional with defaults shown
    'template_prioritization': 'Relevance',
    'precursor_prioritization': 'RelevanceHeuristic',
    'mincount': 0,
    'num_templates': 100,
    'max_cum_prob': 0.995,
    'apply_fast_filter': True,
    'filter_threshold': 0.75,
    
    # modified for this example
    'num_results': 3 # default is 100
}
resp = requests.get(HOST+'/api/retro/', params=params)
pprint(resp.json())
```

    {u'precursors': [{u'necessary_reagent': u'',
                      u'num_examples': 1223,
                      u'plausibility': 0.998188316822052,
                      u'rank': 1,
                      u'score': -0.005976811431568982,
                      u'smiles': u'CN(C)CCCl.OC(c1ccccc1)c1ccccc1',
                      u'smiles_split': [u'CN(C)CCCl', u'OC(c1ccccc1)c1ccccc1'],
                      u'template_score': 0.33462658524513245,
                      u'templates': [u'59c5118c05581eb9f5753c93',
                                     u'59c5118c05581eb9f5753c9d']},
                     {u'necessary_reagent': u'',
                      u'num_examples': 915,
                      u'plausibility': 0.9775225520133972,
                      u'rank': 2,
                      u'score': -0.011262814101429644,
                      u'smiles': u'CN(C)CCO.OC(c1ccccc1)c1ccccc1',
                      u'smiles_split': [u'CN(C)CCO', u'OC(c1ccccc1)c1ccccc1'],
                      u'template_score': 0.1775755137205124,
                      u'templates': [u'59c5118d05581eb9f5753db4',
                                     u'59c511de05581eb9f5758a9b',
                                     u'59c5122205581eb9f575c32d']},
                     {u'necessary_reagent': u'',
                      u'num_examples': 230,
                      u'plausibility': 0.9893513321876526,
                      u'rank': 3,
                      u'score': -0.013463378679132122,
                      u'smiles': u'CN(C)CCO.ClC(c1ccccc1)c1ccccc1',
                      u'smiles_split': [u'CN(C)CCO', u'ClC(c1ccccc1)c1ccccc1'],
                      u'template_score': 0.1485511213541031,
                      u'templates': [u'59c5118e05581eb9f5753df3']}],
     u'request': {u'apply_fast_filter': [u'True'],
                  u'filter_threshold': [u'0.75'],
                  u'max_cum_prob': [u'0.995'],
                  u'mincount': [u'0'],
                  u'num_results': [u'3'],
                  u'num_templates': [u'100'],
                  u'precursor_prioritization': [u'RelevanceHeuristic'],
                  u'target': [u'CN(C)CCOC(c1ccccc1)c1ccccc1'],
                  u'template_prioritization': [u'Relevance']}}


## /api/context/
Given `reactants` and `products`, suggest reaction contexts (catalyst, reagents, solvents, and temperature). The response will have `contexts`, a list of each suggested context in order of recommendation. The maximum number of results that can be returned is 18.


```python
params = {
    'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1', #required
    'products': 'CN(C)CCOC(c1ccccc1)c1ccccc1', #required
    
    'num_results': 5 # default is 10
}
resp = requests.get(HOST+'/api/context/', params=params)
pprint(resp.json())
```

    {u'contexts': [{u'catalyst': u'',
                    u'reagent': u'Cc1ccccc1.[H][N-][H].[Na+]',
                    u'solvent': u'',
                    u'temperature': 94.47986602783203},
                   {u'catalyst': u'',
                    u'reagent': u'c1ccccc1.[H][N-][H].[Na+]',
                    u'solvent': u'',
                    u'temperature': 101.66519927978516},
                   {u'catalyst': u'',
                    u'reagent': u'Cc1ccccc1C.[H][N-][H].[Na+]',
                    u'solvent': u'',
                    u'temperature': 124.40973663330078},
                   {u'catalyst': u'',
                    u'reagent': u'',
                    u'solvent': u'Cc1ccccc1',
                    u'temperature': 109.34732818603516},
                   {u'catalyst': u'',
                    u'reagent': u'',
                    u'solvent': u'c1ccccc1',
                    u'temperature': 102.02490997314453}],
     u'request': {u'num_results': [u'5'],
                  u'products': [u'CN(C)CCOC(c1ccccc1)c1ccccc1'],
                  u'reactants': [u'CN(C)CCCl.OC(c1ccccc1)c1ccccc1']}}


## /api/forward/
Given `reactants`, and optionally `reagents` and `solvent`, suggest probably products. The response will have `outcomes`, with `smiles` of the product, `prob` (probability) of this being the major product, and relative `score` for reach suggestion.


```python
params = {
    'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1', #required
    
    # optional with defaults shown
    'reagents': '',
    'solvent': '',
    
    'num_results': 5 # default is 100
}
resp = requests.get(HOST+'/api/forward/?', params=params)
pprint(resp.json())
```

    {u'outcomes': [{u'prob': 0.9115045729462279,
                    u'rank': 1,
                    u'score': -63.30739974975586,
                    u'smiles': u'CN(C)CCOC(c1ccccc1)c1ccccc1'},
                   {u'prob': 0.08476252957902111,
                    u'rank': 2,
                    u'score': -65.68751525878906,
                    u'smiles': u'c1ccc(Cc2ccccc2)cc1'},
                   {u'prob': 0.001807120522928764,
                    u'rank': 3,
                    u'score': -69.53076171875,
                    u'smiles': u'CN(C)CCC(c1ccccc1)c1ccccc1'},
                   {u'prob': 0.0007843034232317958,
                    u'rank': 4,
                    u'score': -70.3654556274414,
                    u'smiles': u'CN(C)CCC(O)(c1ccccc1)c1ccccc1'},
                   {u'prob': 0.00048256904206722353,
                    u'rank': 5,
                    u'score': -70.85112762451172,
                    u'smiles': u'CCN(C)C'}],
     u'request': {u'num_results': [u'5'],
                  u'reactants': [u'CN(C)CCCl.OC(c1ccccc1)c1ccccc1'],
                  u'reagents': [u''],
                  u'solvent': [u'']}}


## /api/treebuilder/
Given the `smiles` of a target molecule, and various optional settings, resolve possible synthetic pathways terminating according to the stopping criteria provided through settings. The results are structued as a nested dictionary, where the children of chemicals are reactions and the children of reactions are chemicals.


```python
params = {
    'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1', # required
    
    # optional with defaults shown
    'max_depth': 4,
    'max_branching': 25,
    'expansion_time': 60,
    'max_ppg': 10,
    'template_count': 100,
    'max_cum_prob': 0.995,
    'chemical_property_logic': 'none',
    'max_chemprop_c': 0,
    'max_chemprop_n': 0,
    'max_chemprop_o': 0,
    'max_chemprop_h': 0,
    'chemical_popularity_logic': 'none',
    'min_chempop_reactants': 5,
    'min_chempop_products': 5,
    'filter_threshold': 0.75,
    
    'return_first': 'true' # default is false
}
resp = requests.get(HOST+'/api/treebuilder/', params=params)
pprint(resp.json())
```

    {u'request': {u'chemical_popularity_logic': [u'none'],
                  u'chemical_property_logic': [u'none'],
                  u'expansion_time': [u'60'],
                  u'filter_threshold': [u'0.75'],
                  u'max_branching': [u'25'],
                  u'max_chemprop_c': [u'0'],
                  u'max_chemprop_h': [u'0'],
                  u'max_chemprop_n': [u'0'],
                  u'max_chemprop_o': [u'0'],
                  u'max_cum_prob': [u'0.995'],
                  u'max_depth': [u'4'],
                  u'max_ppg': [u'10'],
                  u'min_chempop_products': [u'5'],
                  u'min_chempop_reactants': [u'5'],
                  u'return_first': [u'true'],
                  u'smiles': [u'CN(C)CCOC(c1ccccc1)c1ccccc1'],
                  u'template_count': [u'100']},
     u'trees': [{u'as_product': 44,
                 u'as_reactant': 59,
                 u'children': [{u'children': [{u'as_product': 102,
                                               u'as_reactant': 2438,
                                               u'children': [],
                                               u'id': 1,
                                               u'is_chemical': True,
                                               u'ppg': 1.0,
                                               u'smiles': u'BrC(c1ccccc1)c1ccccc1'},
                                              {u'as_product': 383,
                                               u'as_reactant': 3643,
                                               u'children': [],
                                               u'id': 2,
                                               u'is_chemical': True,
                                               u'ppg': 1.0,
                                               u'smiles': u'CN(C)CCO'}],
                                u'id': 3,
                                u'is_reaction': True,
                                u'necessary_reagent': u'',
                                u'num_examples': 185,
                                u'plausibility': 0.9988686442375183,
                                u'smiles': u'BrC(c1ccccc1)c1ccccc1.CN(C)CCO>>CN(C)CCOC(c1ccccc1)c1ccccc1',
                                u'template_score': 0.07315371185541153,
                                u'tforms': [u'59c5118d05581eb9f5753dbf']}],
                 u'id': 4,
                 u'is_chemical': True,
                 u'ppg': 1.0,
                 u'smiles': u'CN(C)CCOC(c1ccccc1)c1ccccc1'}]}


## /api/template/
Given a template `id` (these are returned with `/api/retro/` precursors) look up template details, such as `reaction_smarts`, and `references` (Reaxys IDs).


```python
params = {
    'id': '59c5300605581eb9f584df3d' # required
}
resp = requests.get(HOST+'/api/template/', params=params)
pprint(resp.json())
```

    {u'request': {u'id': [u'59c5300605581eb9f584df3d']},
     u'template': {u'_id': u'59c5300605581eb9f584df3d',
                   u'chiral': False,
                   u'count': 13,
                   u'dimer_only': False,
                   u'explicit_H': False,
                   u'incompatible_groups': [],
                   u'intra_only': True,
                   u'necessary_reagent': u'',
                   u'reaction_smarts': u'[#7:5]-[C:4](=[O;D1;H0:6])-[c:3]:[c;H0;D3;+0:1](:[#7;a:2])-[N;H0;D3;+0:9](-[C:10])-[c:8]:[#7;a:7]>>Cl-[c;H0;D3;+0:1](:[#7;a:2]):[c:3]-[C:4](-[#7:5])=[O;D1;H0:6].[#7;a:7]:[c:8]-[NH;D2;+0:9]-[C:10]',
                   u'references': [u'4544223',
                                   u'5028471',
                                   u'5028471',
                                   u'5028471',
                                   u'5028471',
                                   u'5029289',
                                   u'8592467',
                                   u'8593425',
                                   u'8593425',
                                   u'23062042',
                                   u'23062042',
                                   u'24072327',
                                   u'38773479']}}


## /api/fast-filter/
Given `reactants` and `products`, return coarse-filter `score` which can be interpreted as a likelihood the reaction may be successful. The resulting score should not be assumed to be linearly proportional to likelihood, but instead should be used to identify bad suggestions with low scores.


```python
params = {
    'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1', # required
    'products': 'CN(C)CCOC(c1ccccc1)c1ccccc1' # required
}
resp = requests.get(HOST+'/api/fast-filter/', params=params)
pprint(resp.json())
```

    {u'request': {u'products': [u'CN(C)CCOC(c1ccccc1)c1ccccc1'],
                  u'reactants': [u'CN(C)CCCl.OC(c1ccccc1)c1ccccc1']},
     u'score': 0.998188316822052}


## /api/scscore/
Given a `smiles` string of a molecule, return the Synthetic Complexity `score`.


```python
params = {
    'smiles': 'OC(c1ccccc1)c1ccccc1' # required
}
resp = requests.get(HOST+'/api/scscore/', params=params)
pprint(resp.json())
```

    {u'request': {u'smiles': [u'OC(c1ccccc1)c1ccccc1']},
     u'score': 1.5127569539402126}


## /api/price/
Given a `smiles` string of a molecule, return the `price` (resolved to integer price per gram values) in the buyables database. A price of 0.0 means the cheical is not in the buyables database.


```python
params = {
    'smiles': 'OC(c1ccccc1)c1ccccc1' # required
}
resp = requests.get(HOST+'/api/price/', params=params)
pprint(resp.json())
```

    {u'price': 1.0, u'request': {u'smiles': [u'OC(c1ccccc1)c1ccccc1']}}


## /api/celery/
Query the status of celery workers on the server. For each queue, the `active` and `available` will sum to the total number of spawned celery workers, where the active workers represent the current number of workers able to complete tasks.


```python
resp = requests.get(HOST+'/api/celery/')
pprint(resp.json())
```

    {u'queues': {u'cr_coordinator': {u'active': 0, u'available': 2},
                 u'cr_network_worker': {u'active': 0, u'available': 2},
                 u'ft_worker': {u'active': 0, u'available': 2},
                 u'sc_coordinator': {u'active': 0, u'available': 2},
                 u'tb_c_worker': {u'active': 7, u'available': 3},
                 u'tb_coordinator_mcts': {u'active': 0, u'available': 2},
                 u'te_coordinator': {u'active': 0, u'available': 2}}}

