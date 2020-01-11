from django.http import JsonResponse
from django.shortcuts import render
from ..globals import scscorer
from ..utils import ajax_error_wrapper

@ajax_error_wrapper
def ajax_scscore_smiles(request):
    print('Got scscore request')
    data = {'err': False}
    smiles = request.GET.get('smiles', None)
    scscore = scscorer.get_score_from_smiles(smiles, noprice=True)
    print(scscore)
    data['html'] = 'Perceived synthetic complexity [1-5]: {:.4f}'.format(scscore)
    return JsonResponse(data)


def scscoring(request):
    return render(request, 'scscoring.html', {})