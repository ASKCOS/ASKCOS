from django.http import JsonResponse
from rdkit import Chem
from celery.exceptions import TimeoutError
from makeit.utilities.contexts import clean_context
from askcos_site.askcos_celery.treeevaluator.scoring_coordinator import evaluate

TIMEOUT = 30

def template_free(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    reactants = request.GET.get('reactants')
    solvent = request.GET.get('solvent', '')
    reagents = request.GET.get('reagents', '')
    num_results = int(request.GET.get('num_results', 100))
    contexts=[clean_context((None, solvent, reagents, '', -1, -1))]

    if not reactants:
        resp['error'] = 'Required parameter "reactants" missing'
        return JsonResponse(resp, status=400)

    rmol = Chem.MolFromSmiles(reactants)
    if not rmol:
        resp['error'] = 'Cannot parse reactants smiles with rdkit'
        return JsonResponse(resp, status=400)

    smol = Chem.MolFromSmiles(contexts[0][1])
    if not smol:
        resp['error'] = 'Cannot parse solvent smiles with rdkit'
        return JsonResponse(resp, status=400)

    remol = Chem.MolFromSmiles(contexts[0][2])
    if not remol:
        resp['error'] = 'Cannot parse reagents smiles with rdkit'
        return JsonResponse(resp, status=400)

    res = evaluate.delay(
        reactants, '', contexts=contexts, forward_scorer='Template_Free', 
        top_n=num_results, return_all_outcomes=True
    )
    
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
    
    outcomes = outcomes[0]['outcomes']
    for out in outcomes:
        o = out.pop('outcome')
        out['smiles'] = o['smiles']
    resp['outcomes'] = outcomes
    return JsonResponse(resp)
