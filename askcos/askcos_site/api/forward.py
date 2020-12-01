from django.http import JsonResponse
from makeit.utilities.contexts import clean_context
from askcos_site.askcos_celery.treeevaluator.scoring_coordinator import evaluate

def template_free(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    reactants = request.GET.get('reactants')
    solvent = request.GET.get('solvent', '')
    reagents = request.GET.get('reagents', '')
    num_results = int(request.GET.get('num_results', 100))
    contexts=[clean_context((None, solvent, reagents, '', -1, -1))]
    res = evaluate.delay(
        reactants, '', contexts=contexts, forward_scorer='Template_Free', 
        top_n=num_results, return_all_outcomes=True
    )
    
    try:
        outcomes = res.get(60)
    except:
        resp['error'] = 'API request timed out (limit 60s)'
        res.revoke()
        return JsonResponse(resp)
    
    outcomes = outcomes[0]['outcomes']
    for out in outcomes:
        o = out.pop('outcome')
        out['smiles'] = o['smiles']
    resp['outcomes'] = outcomes
    return JsonResponse(resp)