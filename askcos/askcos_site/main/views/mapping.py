from django.shortcuts import render, HttpResponse, redirect
from django.urls import reverse
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.conf import settings
import django.contrib.auth.views

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from ..utils import ajax_error_wrapper, resolve_smiles
from ..forms import DrawingInputForm
import time
from askcos_site.askcos_celery.atom_mapper.atom_mapping_worker import get_atom_mapping
import re


# @login_required
def atom_mapping(request):
    return render(request, 'mapping.html')

@ajax_error_wrapper
def ajax_find_atom_mapping(request):
    '''Perform the forward synthesis'''
    data = {'err': False}

    rxnsmiles = request.GET.get('rxnsmiles', '')

    mapper = request.GET.get('mapper', 'WLN atom mapper')

    print('reactants: {}'.format(rxnsmiles))
    print('mapper: {}'.format(mapper))

    startTime = time.time()

    ## need to work on it
    result = get_atom_mapping.delay(rxnsmiles, mapper=mapper)
    data['rxnsmiles_mapped'] = result.get(10)
    print('---------------------------------')
    print(data)
    print('---------------------------------')

    data['html_time'] = '{:.3f} seconds elapsed'.format(time.time() - startTime)

    if data['rxnsmiles_mapped']:
        url = reverse('draw_mapped_reaction', kwargs={'smiles': data['rxnsmiles_mapped']})
        data['html'] = '<img src="' + url + '">'
    else:
        data['html'] = 'No atom mapping is found'

    return JsonResponse(data)
