from django.shortcuts import render, HttpResponse, redirect
from django.template.loader import render_to_string
from django.urls import reverse
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.conf import settings
import django.contrib.auth.views
from pymongo.message import bson
from bson.objectid import ObjectId
import time
import numpy as np
import json
import os
from celery.result import AsyncResult

# TODO: fix this Celery reference
from ...askcos_celery.contextrecommender.cr_coordinator import get_context_recommendations
from askcos_site.askcos_celery.impurity.impurity_worker import get_impurities

from ..globals import PREDICTOR_FOOTNOTE, solvent_choices
from ..utils import ajax_error_wrapper, fix_rgt_cat_slvt, \
    trim_trailing_period
from makeit.utilities.contexts import clean_context
from celery.result import AsyncResult, allow_join_result
from abc import ABCMeta, abstractmethod

#@login_required
def impurity_interactive(request,
                         reactants='', products='', reagents='', solvents='',
                         predictor='WLN forward predictor',
                         inspector='Reaxys inspector',
                         mapper='WLN atom mapper',
                         top_k=3, threshold=0.75):
    '''Builds an impurity page'''

    context = {}
    context['footnote'] = PREDICTOR_FOOTNOTE

    context['reactants'] = reactants
    context['products'] = products
    context['reagents'] = reagents
    context['solvents'] = solvents

    context['predictor'] = predictor
    context['inspector'] = inspector
    context['mapper'] = mapper

    context['top_k'] = str(top_k)
    context['threshold'] = str(threshold)

    return render(request, 'impurity_interactive.html', context)


@ajax_error_wrapper
def ajax_start_impurity(request):
    '''Perform the forward synthesis'''
    data = {'err': False}

    reactants = request.GET.get('reactants', '')
    products = request.GET.get('products', '')
    reagents = request.GET.get('reagents', '')
    solvents = request.GET.get('solvents', '')

    predictor = request.GET.get('predictor', 'WLN forward predictor')
    inspector = request.GET.get('inspector', 'Reaxys inspector')
    mapper = request.GET.get('mapper', 'WLN atom mapper')

    top_k = int(request.GET.get('top_k', 3))
    threshold = float(request.GET.get('threshold', 0.75))

    check_mapping = False
    if request.GET.get('check_mapping', 'True') == 'True':
        check_mapping = True

    print('Conditions for forward synthesis:')
    print('reactants: {}'.format(reactants))
    print('Major products: {}'.format(products))
    print('solvents: {}'.format(solvents))
    print('reagents: {}'.format(reagents))

    print('predictor: {}'.format(predictor))
    print('inspector: {}'.format(inspector))
    print('mapper: {}'.format(mapper))

    print('top_k: {}'.format(top_k))
    print('threshold: {}'.format(threshold))

    print('check mapping {}'.format(check_mapping))


    startTime = time.time()

    ## need to work on it
    result = get_impurities.delay(reactants, reagents=reagents, products=products, solvents=solvents,
                                  predictor_selection=predictor,
                                  inspector_selection=inspector,
                                  mapper_selection=mapper,
                                  top_k=top_k, threshold=threshold, check_mapping=check_mapping)
    data['task_id'] = result.task_id
    print('---------------------------------')
    print(data)
    print('---------------------------------')
    return JsonResponse(data)

@ajax_error_wrapper
def ajax_get_progress(request, task_id):
    output = {}
    print('---------------------------------')
    print(task_id)
    print('---------------------------------')
    result = AsyncResult(task_id)
    state = result.state
    info = result.info
    print(info)
    if state == 'running':
        output['complete'] = False
    else:
        output['complete'] = True
        outcomes = result.get(task_id)
        output['results'] = outcomes
        output['html'] = render_to_string('impurity_outcomes_only.html',
                                        {'outcomes': outcomes['predict_expand']})
        print(outcomes['predict_expand'])
    try:
        output['percent'] = info['percent']
        output['message'] = info['message']
    except :
        output['percent'] = 1
        output['message'] = 'All prediction done!'

    print(output['message'], ' percent: ', output['percent'])
    return HttpResponse(json.dumps(output), content_type='application/json')

