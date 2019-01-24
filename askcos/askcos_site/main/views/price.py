from django.shortcuts import render, HttpResponse, redirect
from django.template.loader import render_to_string
from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.conf import settings
import json

from ..globals import Pricer
from ..utils import ajax_error_wrapper

def price_smiles_func(smiles):
    return Pricer.lookup_smiles(smiles, alreadyCanonical=True)

@login_required
def price_smiles(request, smiles):
    response = HttpResponse(content_type = 'text/plain')
    ppg = Pricer.lookup_smiles(smiles, alreadyCanonical=True)
    if ppg:
        response.write('${}/g'.format(ppg))
    else:
        response.write('cannot buy')
    return response

@ajax_error_wrapper
def ajax_price_smiles(request):
    print('Got price request')
    smiles = request.GET.get('smiles', None)
    data = {'err': False, 'smiles': smiles}
    isomericSmiles = json.loads(request.GET.get('isomericSmiles', 'false'))
    print('isomericSmiles: {}'.format(isomericSmiles))
    data['ppg'] = Pricer.lookup_smiles(smiles, alreadyCanonical=False, isomericSmiles=isomericSmiles)
    print('Result: {}'.format(data['ppg']))
    data['buyable'] = data['ppg'] != 0.
    if data['ppg'] == 0.:
        data['html'] = 'This chemical is <b>not</b> in our database currently'
    else:
        data['html'] = 'This chemical is purchaseable for an estimated <b>$%i/g</b>' % data['ppg']
    return JsonResponse(data)

@login_required 
def pricing(request):
    return render(request, 'pricing.html', {})

@login_required
def price_xrn(request, xrn):
    response = HttpResponse(content_type = 'text/plain')
    ppg = Pricer.lookup_xrn(xrn)
    if ppg:
        response.write('${}/g'.format(ppg))
    else:
        response.write('cannot buy')
    return response