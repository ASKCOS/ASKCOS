from django.shortcuts import render, HttpResponse, redirect
from django.template.loader import render_to_string
from django.urls import reverse
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.conf import settings
import json

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

#@login_required
def scscoring(request):
    return render(request, 'scscoring.html', {})
