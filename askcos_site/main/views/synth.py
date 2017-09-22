from django.shortcuts import render, HttpResponse, redirect
from django.template.loader import render_to_string
from django.core.urlresolvers import reverse
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

from ...askcos_celery.contextrecommender.worker import get_context_recommendation
from ...askcos_celery.forwardpredictor.coordinator import get_outcomes

from ..globals import PREDICTOR_FOOTNOTE, SOLVENT_DB
from ..utils import ajax_error_wrapper, fix_rgt_cat_slvt, \
    trim_trailing_period

@login_required
def synth_interactive(request, reactants='', reagents='', solvent='toluene', 
        temperature='20', mincount='25', product=None):
    '''Builds an interactive forward synthesis page'''

    context = {} 
    context['footnote'] = PREDICTOR_FOOTNOTE
    solvent_choices = []
    for doc in SOLVENT_DB.find({'_id': {'$ne': 'default'}}):
        solvent_choices.append({
            'smiles': doc['smiles'],
            'name': doc['name'],
        })
    context['solvent_choices'] = sorted(solvent_choices, key = lambda x: x['name'])
    context['reactants'] = reactants

    if product is None:    
        context['reagents'] = reagents
        context['solventselected'] = solvent
        context['temperature'] = temperature
    else:
        # Get suggested conditions
        smiles = '%s>>%s' % (reactants, product)
        res = get_context_recommendation.delay(smiles, n=10)
        contexts = res.get(60)
        if contexts is None or len(contexts) == 0:
            raise ValueError('Context recommender was unable to get valid context(?)')
        (T1, slvt1, rgt1, cat1, t1, y1) = contexts[0]
        slvt1 = trim_trailing_period(slvt1)
        rgt1 = trim_trailing_period(rgt1)
        cat1 = trim_trailing_period(cat1)
        (rgt1, cat1, slvt1) = fix_rgt_cat_slvt(rgt1, cat1, slvt1)
        for slvt in solvent_choices:
            if slvt['smiles'] == slvt1:
                context['solventselected'] = slvt['name']
                break
        context['reagents'] = rgt1
        context['temperature'] = T1 

    context['mincount'] = mincount if mincount != '' else settings.SYNTH_TRANSFORMS['mincount']
    return render(request, 'synth_interactive.html', context)

@login_required
def synth_interactive_smiles(request, smiles):
    '''Synth interactive initialized w/ reaction smiles'''
    return synth_interactive(request, reactants=smiles.split('>')[0], product=smiles.split('>')[-1])

@ajax_error_wrapper
def ajax_validate_temperature_synth(request):
    '''Checks to see if a temperature is valid'''
    data = {'err': False}

    temperature = request.GET.get('temperature', None)
    print('temperature from Ajax: {}'.format(temperature))

    try:
        temperature = float(temperature)
        data['temperature'] = temperature
    except Exception as e:
        data['err'] = True

    return JsonResponse(data)

@ajax_error_wrapper
def ajax_start_synth(request):
    '''Perform the forward synthesis'''
    data = {'err': False}

    reactants = request.GET.get('reactants', None)
    solvent = request.GET.get('solvent', None)
    temperature = request.GET.get('temperature', None)
    reagents = request.GET.get('reagents', None)
    mincount = int(request.GET.get('mincount', None))
    maxreturn = int(request.GET.get('maxreturn', 100))
    print('Conditions for forward synthesis:')
    print('reactants: {}'.format(reactants))
    print('solvent: {}'.format(solvent))
    print('temp: {}'.format(temperature))
    print('reagents: {}'.format(reagents))
    print('mincount: {}'.format(mincount))
    print('max return: {}'.format(maxreturn))

    startTime = time.time()
    from askcos_site.askcos_celery.forwardpredictor.coordinator import get_outcomes
    res = get_outcomes.delay(reactants, 
        contexts=[(temperature, reagents, solvent)], 
        mincount=mincount,
        top_n=maxreturn)
    outcomes = res.get(300)[0]

    print('Got top outcomes, length {}'.format(len(outcomes)))
    data['html_time'] = '{:.3f} seconds elapsed'.format(time.time() - startTime)

    if outcomes:

        data['html'] = render_to_string('synth_outcomes_only.html', 
            {'outcomes': outcomes})
    else:
        data['html'] = 'No outcomes found? That is weird...'
    
    # Save in session in case used wants to print
    request.session['last_synth_interactive'] = {'reactants': reactants, 
        'temperature': temperature, 'reagents': reagents, 'solvent': solvent,
        'mincount': mincount, 'outcomes': outcomes}

    return JsonResponse(data)

