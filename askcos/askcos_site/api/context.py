from django.http import JsonResponse
from celery.exceptions import TimeoutError
from askcos_site.askcos_celery.contextrecommender.cr_nn_worker import get_n_conditions as neighbor_get_n_conditions
from askcos_site.askcos_celery.contextrecommender.cr_network_worker import get_n_conditions as network_get_n_conditions
from rdkit import Chem

TIMEOUT = 30

TIMEOUT = 30

def neural_network(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    reactants = request.GET.get('reactants')
    products = request.GET.get('products')
    with_smiles = request.GET.get('with_smiles', 'True') in ['True', 'true']
    singleSlvt = request.GET.get('singleSlvt', 'True') in ['True', 'true']
    return_scores = request.GET.get('return_scores') in ['True', 'true']
    n = int(request.GET.get('num_results', 10))
    rxn = reactants+'>>'+products

    if not reactants:
        resp['error'] = 'Required parameter "reactants" missing'
        return JsonResponse(resp, status=400)
    
    if not products:
        resp['error'] = 'Required parameter "products" missing'
        return JsonResponse(resp, status=400)

    rmol = Chem.MolFromSmiles(reactants)
    if not rmol:
        resp['error'] = 'Cannot parse reactants smiles with rdkit'
        return JsonResponse(resp, status=400)

    pmol = Chem.MolFromSmiles(products)
    if not pmol:
        resp['error'] = 'Cannot parse products smiles with rdkit'
        return JsonResponse(resp, status=400)

    res = network_get_n_conditions.delay(rxn, n, singleSlvt, with_smiles, return_scores)
    try:
        if return_scores:
            contexts, scores = res.get(TIMEOUT)
        else:
            contexts = res.get(TIMEOUT)
    except TimeoutError:
        resp['error'] = 'API request timed out (limit {}s)'.format(TIMEOUT)
        res.revoke()
        return JsonResponse(resp, status=408)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return JsonResponse(resp, status=400)

    json_contexts = []
    for context in contexts:
        c = {
            'temperature': context[0],
            'solvent': context[1],
            'reagent': context[2],
            'catalyst': context[3]
        }
        json_contexts.append(c)
    if return_scores:
        for c, s in zip(json_contexts, scores):
            c['score'] = s
    resp['contexts'] = json_contexts
    return JsonResponse(resp)
