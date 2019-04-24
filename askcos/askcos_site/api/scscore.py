from django.http import JsonResponse
from makeit.prioritization.precursors.scscore import SCScorePrecursorPrioritizer

scscorer = SCScorePrecursorPrioritizer()
scscorer.load_model(model_tag='1024bool')

def scscore(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    smiles = request.GET.get('smiles', None)
    scscore = scscorer.get_score_from_smiles(smiles, noprice=True)
    resp['score'] = scscore
    return JsonResponse(resp)