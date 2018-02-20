from django.shortcuts import render, HttpResponse, redirect
from django.template.loader import render_to_string
from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.conf import settings
import django.contrib.auth.views
from datetime import datetime
import time
import numpy as np 
import json
import os
import rdkit.Chem as Chem 
import makeit.global_config as makeit_gc

from ...askcos_celery.contextrecommender.cr_coordinator import get_context_recommendations

from ..utils import ajax_error_wrapper, fix_rgt_cat_slvt, \
    trim_trailing_period

@login_required
def context_rxnsmiles(request, smiles=None, reactants=None, product=None):
    context = {}
    if smiles is not None:
        context['reactants'] = smiles.split('>')[0]
        context['product'] = smiles.split('>')[-1]
    else:
        if reactants is not None:
            context['reactants'] = reactants
        if product is not None:
            context['product'] = product
    return render(request, 'context.html', context)

@ajax_error_wrapper
def ajax_context_rxnsmiles(request):
    '''Evaluate rxn_smiles'''
    data = {'err': False}
    smiles = request.GET.get('smiles', None)
    verbose = json.loads(request.GET.get('verbose', 'false'))
    context_recommender = request.GET.get('context_recommender', 'Neural_Network')

    reactants = smiles.split('>>')[0].split('.')
    products = smiles.split('>>')[1].split('.')
    print('...trying to get predicted context')

    res = get_context_recommendations.delay(smiles, n=10, singleSlvt=False,
        context_recommender=context_recommender)
    contexts = res.get(60)
    print('Got context(s)')
    print(contexts)
    if contexts is None:
        raise ValueError('Context recommender was unable to get valid context(?)')

    if contexts:
        data['html'] = render_to_string('context_recs_only.html', 
            {'contexts': [context_to_dict(x) for x in contexts],
             'reactants': '.'.join(reactants)})
    else:
        data['html'] = 'No recommendations found? That is weird...'

    return JsonResponse(data)

def context_to_dict(context):
    (T1, slvt1, rgt1, cat1, t1, y1) = context
    return {
        'temperature': T1,
        'solvents': slvt1 if slvt1 != '.' else '',
        'reagents': rgt1,
        'reagents_combined': '.'.join([x for x in rgt1.split('.') + cat1.split('.') if x]),
        'catalysts': cat1,
        'time': t1,
        'yield': y1,
    }

@login_required
def context_rxnsmiles_target(request, smiles):
    '''Synth interactive initialized w/ reaction smiles'''
    return context_rxnsmiles(request, reactants=smiles.split('>')[0], product=smiles.split('>')[-1])

@login_required
def context_rxnsmiles_target2(request, reactants, product):
    '''Synth interactive initialized w/ reaction smiles'''
    return context_rxnsmiles(request, reactants=reactants, product=product)