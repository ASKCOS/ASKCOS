from __future__ import absolute_import, unicode_literals, print_function
from django.shortcuts import render 
from django.http import JsonResponse 
from django.conf import settings
from django.template.loader import render_to_string
from ...askcos_celery.siteselectivity.sites_worker import get_sites
from ..utils import ajax_error_wrapper 

def site_prediction(request, err=None):
    '''
    Site Predictor webpage
    '''
    return render(request, 'sites.html', {'err': err})

@ajax_error_wrapper
def ajax_get_sites(request):
    '''Evaluate rxn_smiles'''
    data = {'err': False}
    smiles = request.GET.get('smiles', None)
    
    res = get_sites.delay(smiles)
    try:
        result = res.get(20)
    except:
        raise ValueError('Timeout for some crazy reason...')
    if res is None:
        raise ValueError('Recommender was unable to get valid result(?)')
    
    if result:
        data['html'] = render_to_string('sites_recs_only.html', 
            {'smiles': smiles,
             'results': result
             })

    else:
        data['html'] = 'No recommendations found? That is weird...'

    return JsonResponse(data)