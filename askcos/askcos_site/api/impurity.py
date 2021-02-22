from django.http import JsonResponse
from rdkit import Chem
from celery.exceptions import TimeoutError
from makeit.utilities.contexts import clean_context
# TODO: need to work on it
from askcos_site.askcos_celery.impurity.impurity_worker import get_outcomes

TIMEOUT = 30


def impurity_predict(request):
    resp = {}
    resp['request'] = dict(**request.GET)

    reactants = request.GET.get('reactants')
    solvents = request.GET.get('solvents', '')
    reagents = request.GET.get('reagents', '')
    products = request.GET.get('prodcuts', '')

    # TODO: add model selections

    if not reactants:
        resp['error'] = 'Required parameter "reactants" missing'
        return JsonResponse(resp, status=400)

    rmol = Chem.MolFromSmiles(reactants)
    if not rmol:
        resp['error'] = 'Cannot parse reactants smiles with rdkit'
        return JsonResponse(resp, status=400)

    smol = Chem.MolFromSmiles(solvents)
    if not smol:
        resp['error'] = 'Cannot parse solvent smiles with rdkit'
        return JsonResponse(resp, status=400)

    remol = Chem.MolFromSmiles(reagents)
    if not remol:
        resp['error'] = 'Cannot parse reagents smiles with rdkit'
        return JsonResponse(resp, status=400)

    pmol = Chem.MolFromSmiles(products)
    if not pmol:
        resp['error'] = 'Cannot parse products smiles with rdkit'
        return JsonResponse(resp, status=400)

    # TODO: need work
    res = get_outcomes.delay(
        reactants, reagents=reagents, products=products, solvents=solvents)

    try:
        outcomes = res.get(TIMEOUT)
    except TimeoutError:
        resp['error'] = 'API request timed out (limit {}s)'.format(TIMEOUT)
        res.revoke()
        return JsonResponse(resp, status=408)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return JsonResponse(resp, status=400)

    # TODO: need work
    outcomes = outcomes[0]['outcomes']
    for out in outcomes:
        o = out.pop('outcome')
        out['smiles'] = o['smiles']
    resp['outcomes'] = outcomes
    return JsonResponse(resp)
