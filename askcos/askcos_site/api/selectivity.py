from rdkit import Chem
from django.http import JsonResponse
from celery.exceptions import TimeoutError
from askcos_site.askcos_celery.siteselectivity.sites_worker import get_sites

def selectivity(request):
    '''Evaluate rxn_smiles'''
    resp = {}
    smiles = request.GET.get('smiles', None)
    if not smiles:
        resp['error'] = 'Required parameter "smiles" missing'
        return JsonResponse(resp, status=400)

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        resp['error'] = 'Cannot parse smiles smiles with rdkit'
        return JsonResponse(resp, status=400)

    res = get_sites.delay(smiles)
    try:
        result = res.get(30)
    except TimeoutError:
        resp['error'] = 'API request timed out (limit 30s)'
        res.revoke()
        return JsonResponse(resp, status=408)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return JsonResponse(resp, status=400)

    #if result:
    resp['results'] = {'smiles': smiles,
            'task_scores': result}

    return JsonResponse(resp)

