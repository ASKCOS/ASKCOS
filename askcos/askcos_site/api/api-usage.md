

```python
import requests
from pprint import pprint
```


```python
HOST = 'https://<your-askcos-ip>'
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

The API now sits at a server that expects communication using HTTPS. By default, a randomly generated self-signed certificate is used, which is considered insecure (despite being much more secure than not using HTTPS at all). The `requests` Python package, by default, will not allow you to use HTTPS to communicate with a server that doesn't have a valid certificate, however this behavior can be overriden by passing `verify=False` to the get request, as shown in the examples below.

Internally, the code behind the API will be converting any SMILES string you provide to a canonicalized SMILES string. Therefore, the results returned by the API may include a slightly modified string than the result parameters you provided. Keep note of this if using SMILES string searching to parse through results.

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
resp = requests.get(HOST+'/api/retro/', params=params, verify=False)
pprint(resp.json())
```

    /usr/local/lib/python3.5/dist-packages/urllib3/connectionpool.py:847: InsecureRequestWarning: Unverified HTTPS request is being made. Adding certificate verification is strongly advised. See: https://urllib3.readthedocs.io/en/latest/advanced-usage.html#ssl-warnings
      InsecureRequestWarning)


    {'precursors': [{'necessary_reagent': '',
                     'num_examples': 1223,
                     'plausibility': 0.998188316822052,
                     'rank': 1,
                     'score': -0.005976811431568982,
                     'smiles': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1',
                     'smiles_split': ['CN(C)CCCl', 'OC(c1ccccc1)c1ccccc1'],
                     'template_score': 0.33462658524513245,
                     'templates': ['59c5118c05581eb9f5753c93',
                                   '59c5118c05581eb9f5753c9d']},
                    {'necessary_reagent': '',
                     'num_examples': 915,
                     'plausibility': 0.9775225520133972,
                     'rank': 2,
                     'score': -0.011262814101429644,
                     'smiles': 'CN(C)CCO.OC(c1ccccc1)c1ccccc1',
                     'smiles_split': ['CN(C)CCO', 'OC(c1ccccc1)c1ccccc1'],
                     'template_score': 0.1775755137205124,
                     'templates': ['59c5118d05581eb9f5753db4',
                                   '59c511de05581eb9f5758a9b',
                                   '59c5122205581eb9f575c32d']},
                    {'necessary_reagent': '',
                     'num_examples': 230,
                     'plausibility': 0.9893513321876526,
                     'rank': 3,
                     'score': -0.013463378679132122,
                     'smiles': 'CN(C)CCO.ClC(c1ccccc1)c1ccccc1',
                     'smiles_split': ['CN(C)CCO', 'ClC(c1ccccc1)c1ccccc1'],
                     'template_score': 0.1485511213541031,
                     'templates': ['59c5118e05581eb9f5753df3']}],
     'request': {'apply_fast_filter': ['True'],
                 'filter_threshold': ['0.75'],
                 'max_cum_prob': ['0.995'],
                 'mincount': ['0'],
                 'num_results': ['3'],
                 'num_templates': ['100'],
                 'precursor_prioritization': ['RelevanceHeuristic'],
                 'target': ['CN(C)CCOC(c1ccccc1)c1ccccc1'],
                 'template_prioritization': ['Relevance']}}


## /api/context/
Given `reactants` and `products`, suggest reaction contexts (catalyst, reagents, solvents, and temperature). The response will have `contexts`, a list of each suggested context in order of recommendation. The maximum number of results that can be returned is 18.


```python
params = {
    'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1', #required
    'products': 'CN(C)CCOC(c1ccccc1)c1ccccc1', #required
    
    'num_results': 5, # default is 10
    'return_scores': 'true' # default is false
}
resp = requests.get(HOST+'/api/context/', params=params, verify=False)
pprint(resp.json())
```

    /usr/local/lib/python3.5/dist-packages/urllib3/connectionpool.py:847: InsecureRequestWarning: Unverified HTTPS request is being made. Adding certificate verification is strongly advised. See: https://urllib3.readthedocs.io/en/latest/advanced-usage.html#ssl-warnings
      InsecureRequestWarning)


    {'contexts': [{'catalyst': '',
                   'reagent': 'Cc1ccccc1.[H][N-][H].[Na+]',
                   'score': 0.3389343023300171,
                   'solvent': '',
                   'temperature': 94.4798812866211},
                  {'catalyst': '',
                   'reagent': 'c1ccccc1.[H][N-][H].[Na+]',
                   'score': 0.12604430317878723,
                   'solvent': '',
                   'temperature': 101.66520690917969},
                  {'catalyst': '',
                   'reagent': 'Cc1ccccc1C.[H][N-][H].[Na+]',
                   'score': 0.10769638419151306,
                   'solvent': '',
                   'temperature': 124.40973663330078},
                  {'catalyst': '',
                   'reagent': '',
                   'score': 0.004865644965320826,
                   'solvent': 'Cc1ccccc1',
                   'temperature': 109.34728240966797},
                  {'catalyst': '',
                   'reagent': '',
                   'score': 0.004308402072638273,
                   'solvent': 'c1ccccc1',
                   'temperature': 102.02490234375}],
     'request': {'num_results': ['5'],
                 'products': ['CN(C)CCOC(c1ccccc1)c1ccccc1'],
                 'reactants': ['CN(C)CCCl.OC(c1ccccc1)c1ccccc1'],
                 'return_scores': ['true']}}


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
resp = requests.get(HOST+'/api/forward/?', params=params, verify=False)
pprint(resp.json())
```

    /usr/local/lib/python3.5/dist-packages/urllib3/connectionpool.py:847: InsecureRequestWarning: Unverified HTTPS request is being made. Adding certificate verification is strongly advised. See: https://urllib3.readthedocs.io/en/latest/advanced-usage.html#ssl-warnings
      InsecureRequestWarning)


    {'outcomes': [{'prob': 0.9115045620575345,
                   'rank': 1,
                   'score': -63.30739974975586,
                   'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1'},
                  {'prob': 0.08476252470830846,
                   'rank': 2,
                   'score': -65.68264284833658,
                   'smiles': 'c1ccc(Cc2ccccc2)cc1'},
                  {'prob': 0.0018071342418599589,
                   'rank': 3,
                   'score': -69.53075408935547,
                   'smiles': 'CN(C)CCC(c1ccccc1)c1ccccc1'},
                  {'prob': 0.0007843034582285129,
                   'rank': 4,
                   'score': -70.3654556274414,
                   'smiles': 'CN(C)CCC(O)(c1ccccc1)c1ccccc1'},
                  {'prob': 0.00048256904120668755,
                   'rank': 5,
                   'score': -70.85112762451172,
                   'smiles': 'CCN(C)C'}],
     'request': {'num_results': ['5'],
                 'reactants': ['CN(C)CCCl.OC(c1ccccc1)c1ccccc1'],
                 'reagents': [''],
                 'solvent': ['']}}


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
resp = requests.get(HOST+'/api/treebuilder/', params=params, verify=False)
pprint(resp.json())
```

    /usr/local/lib/python3.5/dist-packages/urllib3/connectionpool.py:847: InsecureRequestWarning: Unverified HTTPS request is being made. Adding certificate verification is strongly advised. See: https://urllib3.readthedocs.io/en/latest/advanced-usage.html#ssl-warnings
      InsecureRequestWarning)


    {'request': {'chemical_popularity_logic': ['none'],
                 'chemical_property_logic': ['none'],
                 'expansion_time': ['60'],
                 'filter_threshold': ['0.75'],
                 'max_branching': ['25'],
                 'max_chemprop_c': ['0'],
                 'max_chemprop_h': ['0'],
                 'max_chemprop_n': ['0'],
                 'max_chemprop_o': ['0'],
                 'max_cum_prob': ['0.995'],
                 'max_depth': ['4'],
                 'max_ppg': ['10'],
                 'min_chempop_products': ['5'],
                 'min_chempop_reactants': ['5'],
                 'return_first': ['true'],
                 'smiles': ['CN(C)CCOC(c1ccccc1)c1ccccc1'],
                 'template_count': ['100']},
     'trees': [{'as_product': 44,
                'as_reactant': 59,
                'children': [{'children': [{'as_product': 172,
                                            'as_reactant': 4651,
                                            'children': [],
                                            'id': 1,
                                            'is_chemical': True,
                                            'ppg': 1.0,
                                            'smiles': 'CN(C)CCCl'},
                                           {'as_product': 2004,
                                            'as_reactant': 5783,
                                            'children': [],
                                            'id': 2,
                                            'is_chemical': True,
                                            'ppg': 1.0,
                                            'smiles': 'OC(c1ccccc1)c1ccccc1'}],
                              'id': 3,
                              'is_reaction': True,
                              'necessary_reagent': '',
                              'num_examples': 607,
                              'plausibility': 0.998188316822052,
                              'smiles': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1>>CN(C)CCOC(c1ccccc1)c1ccccc1',
                              'template_score': 0.33462631702423096,
                              'tforms': ['59c5118c05581eb9f5753c93']}],
                'id': 4,
                'is_chemical': True,
                'ppg': 1.0,
                'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1'},
               {'as_product': 44,
                'as_reactant': 59,
                'children': [{'children': [{'as_product': 383,
                                            'as_reactant': 3643,
                                            'children': [],
                                            'id': 5,
                                            'is_chemical': True,
                                            'ppg': 1.0,
                                            'smiles': 'CN(C)CCO'},
                                           {'as_product': 2004,
                                            'as_reactant': 5783,
                                            'children': [],
                                            'id': 2,
                                            'is_chemical': True,
                                            'ppg': 1.0,
                                            'smiles': 'OC(c1ccccc1)c1ccccc1'}],
                              'id': 6,
                              'is_reaction': True,
                              'necessary_reagent': '',
                              'num_examples': 266,
                              'plausibility': 0.9775225520133972,
                              'smiles': 'CN(C)CCO.OC(c1ccccc1)c1ccccc1>>CN(C)CCOC(c1ccccc1)c1ccccc1',
                              'template_score': 0.177575021982193,
                              'tforms': ['59c5118d05581eb9f5753db4']}],
                'id': 4,
                'is_chemical': True,
                'ppg': 1.0,
                'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1'}]}


## /api/template/
Given a template `id` (these are returned with `/api/retro/` precursors) look up template details, such as `reaction_smarts`, and `references` (Reaxys IDs).


```python
params = {
    'id': '59c5300605581eb9f584df3d' # required
}
resp = requests.get(HOST+'/api/template/', params=params, verify=False)
pprint(resp.json())
```

    {'request': {'id': ['59c5300605581eb9f584df3d']},
     'template': {'_id': '59c5300605581eb9f584df3d',
                  'count': 13,
                  'dimer_only': False,
                  'intra_only': True,
                  'necessary_reagent': '',
                  'reaction_smarts': '[#7:5]-[C:4](=[O;D1;H0:6])-[c:3]:[c;H0;D3;+0:1](:[#7;a:2])-[N;H0;D3;+0:9](-[C:10])-[c:8]:[#7;a:7]>>Cl-[c;H0;D3;+0:1](:[#7;a:2]):[c:3]-[C:4](-[#7:5])=[O;D1;H0:6].[#7;a:7]:[c:8]-[NH;D2;+0:9]-[C:10]',
                  'references': ['4544223',
                                 '5028471',
                                 '5028471',
                                 '5028471',
                                 '5028471',
                                 '5029289',
                                 '8592467',
                                 '8593425',
                                 '8593425',
                                 '23062042',
                                 '23062042',
                                 '24072327',
                                 '38773479']}}


    /usr/local/lib/python3.5/dist-packages/urllib3/connectionpool.py:847: InsecureRequestWarning: Unverified HTTPS request is being made. Adding certificate verification is strongly advised. See: https://urllib3.readthedocs.io/en/latest/advanced-usage.html#ssl-warnings
      InsecureRequestWarning)


## /api/fast-filter/
Given `reactants` and `products`, return coarse-filter `score` which can be interpreted as a likelihood the reaction may be successful. The resulting score should not be assumed to be linearly proportional to likelihood, but instead should be used to identify bad suggestions with low scores.


```python
params = {
    'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1', # required
    'products': 'CN(C)CCOC(c1ccccc1)c1ccccc1' # required
}
resp = requests.get(HOST+'/api/fast-filter/', params=params, verify=False)
pprint(resp.json())
```

    /usr/local/lib/python3.5/dist-packages/urllib3/connectionpool.py:847: InsecureRequestWarning: Unverified HTTPS request is being made. Adding certificate verification is strongly advised. See: https://urllib3.readthedocs.io/en/latest/advanced-usage.html#ssl-warnings
      InsecureRequestWarning)


    {'request': {'products': ['CN(C)CCOC(c1ccccc1)c1ccccc1'],
                 'reactants': ['CN(C)CCCl.OC(c1ccccc1)c1ccccc1']},
     'score': 0.998188316822052}


## /api/scscore/
Given a `smiles` string of a molecule, return the Synthetic Complexity `score`.


```python
params = {
    'smiles': 'OC(c1ccccc1)c1ccccc1' # required
}
resp = requests.get(HOST+'/api/scscore/', params=params, verify=False)
pprint(resp.json())
```

    {'request': {'smiles': ['OC(c1ccccc1)c1ccccc1']}, 'score': 1.5127569539402126}


    /usr/local/lib/python3.5/dist-packages/urllib3/connectionpool.py:847: InsecureRequestWarning: Unverified HTTPS request is being made. Adding certificate verification is strongly advised. See: https://urllib3.readthedocs.io/en/latest/advanced-usage.html#ssl-warnings
      InsecureRequestWarning)


## /api/price/
Given a `smiles` string of a molecule, return the `price` (resolved to integer price per gram values) in the buyables database. A price of 0.0 means the cheical is not in the buyables database.


```python
params = {
    'smiles': 'OC(c1ccccc1)c1ccccc1' # required
}
resp = requests.get(HOST+'/api/price/', params=params, verify=False)
pprint(resp.json())
```

    /usr/local/lib/python3.5/dist-packages/urllib3/connectionpool.py:847: InsecureRequestWarning: Unverified HTTPS request is being made. Adding certificate verification is strongly advised. See: https://urllib3.readthedocs.io/en/latest/advanced-usage.html#ssl-warnings
      InsecureRequestWarning)


    {'price': 1.0, 'request': {'smiles': ['OC(c1ccccc1)c1ccccc1']}}


## /api/celery/
Query the status of celery workers on the server. For each queue, the `active` and `available` will sum to the total number of spawned celery workers, where the active workers represent the current number of workers able to complete tasks.


```python
resp = requests.get(HOST+'/api/celery/', verify=False)
pprint(resp.json())
```

    /usr/local/lib/python3.5/dist-packages/urllib3/connectionpool.py:847: InsecureRequestWarning: Unverified HTTPS request is being made. Adding certificate verification is strongly advised. See: https://urllib3.readthedocs.io/en/latest/advanced-usage.html#ssl-warnings
      InsecureRequestWarning)


    {'queues': [{'available': 2,
                 'busy': 0,
                 'name': 'Context Recommender Coordinator',
                 'queue': 'cr_coordinator'},
                {'available': 2,
                 'busy': 0,
                 'name': 'Context Recommender Worker',
                 'queue': 'cr_network_worker'},
                {'available': 2,
                 'busy': 0,
                 'name': 'Forward Predictor Scoring Coordinator',
                 'queue': 'sc_coordinator'},
                {'available': 2,
                 'busy': 0,
                 'name': 'Forward Predictor Worker',
                 'queue': 'ft_worker'},
                {'available': 10,
                 'busy': 0,
                 'name': 'One-Step/Tree Builder Retrosynthesis Worker',
                 'queue': 'tb_c_worker'},
                {'available': 2,
                 'busy': 0,
                 'name': 'Tree Builder Coordinator',
                 'queue': 'tb_coordinator_mcts'},
                {'available': 2,
                 'busy': 0,
                 'name': 'Tree Evaluation Coordinator',
                 'queue': 'te_coordinator'}]}

