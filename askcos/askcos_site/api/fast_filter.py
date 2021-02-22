from django.http import JsonResponse
from rdkit import Chem
from makeit import global_config as gc
from askcos_site.askcos_celery.treebuilder.tb_c_worker import fast_filter_check

TIMEOUT = 30

def fast_filter(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    reactants = request.GET.get('reactants')
    products = request.GET.get('products')

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

    res = fast_filter_check.delay(reactants, products)
    try:
        outcome = res.get(TIMEOUT)
    except TimeoutError:
        resp['error'] = 'API request timed out (limit {}s)'.format(TIMEOUT)
        res.revoke()
        return JsonResponse(resp, status=408)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return JsonResponse(resp, status=400)
    
    resp['score'] = outcome
    return JsonResponse(resp)
