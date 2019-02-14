from django.http import JsonResponse
from makeit.utilities.contexts import clean_context
from askcos_site.askcos_celery.treeevaluator.scoring_coordinator import evaluate

def template_free(request):
    resp = {}
    reactants = request.GET.get('reactants')
    solvent = request.GET.get('solvent', '')
    temperature = request.GET.get('temperature', None)
    reagents = request.GET.get('reagents', '')
    maxreturn = int(request.GET.get('maxreturn', 100))
    contexts=[clean_context((temperature, solvent, reagents, '', -1, -1))]
    print(reagents)
    res = evaluate.delay(reactants, '',
        contexts=contexts, 
        forward_scorer='Template_Free', top_n=maxreturn, return_all_outcomes=True)
    outcomes = res.get(300)[0]['outcomes']
    resp['products'] = [o['outcome'] for o in outcomes]
    return JsonResponse(resp)