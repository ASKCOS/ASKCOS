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

# TODO: fix this Celery reference
from ...askcos_celery.contextrecommender.cr_coordinator import get_context_recommendations
from ...askcos_celery.treeevaluator.scoring_coordinator import evaluate

from ..globals import DONE_SYNTH_PREDICTIONS
from ..utils import ajax_error_wrapper, fix_rgt_cat_slvt, \
    trim_trailing_period

@login_required
def evaluate_rxnsmiles(request):
    return render(request, 'evaluate.html', {})

@ajax_error_wrapper
def ajax_evaluate_rxnsmiles(request):
    '''Evaluate rxn_smiles'''
    data = {'err': False}
    smiles = request.GET.get('smiles', None)
    verbose = json.loads(request.GET.get('verbose', 'false'))
    synth_mincount = int(request.GET.get('synth_mincount', 0))
    necessary_reagent = request.GET.get('necessary_reagent', '')
    if necessary_reagent == 'false':
        necessary_reagent = ''
    if '{}-{}'.format(smiles, synth_mincount) in DONE_SYNTH_PREDICTIONS:
        data = DONE_SYNTH_PREDICTIONS['{}-{}'.format(smiles, synth_mincount)]
        return JsonResponse(data)
    reactants = smiles.split('>>')[0].split('.')
    products = smiles.split('>>')[1].split('.')
    print('...trying to get predicted context')

    if necessary_reagent:
        num_contexts = 1
    else:
        num_contexts = 10
    res = get_context_recommendations.delay(smiles, n=num_contexts, context_recommender='Nearest_Neighbor')
    contexts = res.get(60)
    print('Got context(s)')
    print(contexts)
    if contexts is None:
        raise ValueError('Context recommender was unable to get valid context(?)')

    # Run
    reactant_smiles = smiles.split('>>')[0]
    print('Running forward evaluator on {}'.format(reactant_smiles))
    if necessary_reagent and contexts[0][2]:
        reactant_smiles += contexts[0][2] # add rgt
    
    res = evaluate.delay(reactant_smiles, products[0], contexts, 
        forward_scorer='Template_Based', mincount=synth_mincount, top_n=50)
    all_outcomes = res.get(300)


    if all([len(outcome) == 0 for outcome in all_outcomes]):
        if not verbose:
            data['html'] = 'Could not get outcomes - recommended context(s) unparseable'
            for i, (T, slvt, rgt, cat, t, y) in enumerate(contexts):
                data['html'] += '<br>{}) T={}, rgt={}, slvt={}'.format(i+1, T, rgt, slvt)
            data['html_color'] = str('#%02x%02x%02x' % (int(255), int(0), int(0)))
            return JsonResponse(data)
        else:
            # TODO: expand
            data['html'] = '<h3>Could not get outcomes - recommended context(s) unparseable</h3>\n<ol>\n'
            for i, (T, slvt, rgt, cat, t, y) in enumerate(contexts):
                data['html'] += '<li>Temp: {} C<br>Reagents: {}<br>Solvent: {}</li>\n'.format(T, rgt, slvt)
            data['html'] += '</ol>'
            data['html_color'] = str('#%02x%02x%02x' % (int(255), int(0), int(0)))
            return JsonResponse(data)
    plausible = [outcome['target']['prob'] for outcome in all_outcomes]
    ranks = [outcome['target']['rank'] for outcome in all_outcomes]
    major_prods = [outcome['top_product']['smiles'] for outcome in all_outcomes]
    major_probs = [outcome['top_product']['prob'] for outcome in all_outcomes]
    
    best_context_i = np.argmax(plausible)
    plausible = plausible[best_context_i]
    rank = ranks[best_context_i]
    best_context = contexts[best_context_i]
    major_prod = major_prods[best_context_i]
    major_prob = major_probs[best_context_i]

    # Report
    print('Recommended context(s): {}'.format(best_context))
    print('Plausibility: {}'.format(plausible))
    (T1, slvt1, rgt1, cat1, t1, y1) = best_context

    if not verbose:
        if not rgt1: rgt1 = 'no '
        data['html'] = 'Plausibility score: {} (rank {})'.format(plausible, rank)
        data['html'] += '<br><br><u>Top conditions</u>'
        data['html'] += '<br>{} C'.format(T1)
        data['html'] += '<br>{} solvent'.format(slvt1)
        data['html'] += '<br>{} reagents'.format(rgt1)
        data['html'] += '<br>nearest-neighbor got {}% yield'.format(y1)
        if rank != 1:
            data['html'] += '<br>Predicted major product with p = {}'.format(major_prob)
            data['html'] += '<br>{}'.format(major_prod)
            if major_prod != 'none found':
                url = reverse('draw_smiles', kwargs={'smiles':major_prod})
                data['html'] += '<br><img src="' + url + '">'
        data['html'] += '<br>(calc. used synth_mincount {})'.format(synth_mincount)
    else:
        if not rgt1: rgt1 = 'none'
        data['html'] = '<h3>Plausibility score: {} (rank {})</h3>'.format(plausible, rank)
        data['html'] += '\n<br><u>Proposed conditions ({} tried)</u>\n'.format(len(contexts_for_predictor))
        data['html'] += '<br>Temp: {} C<br>Reagents: {}<br>Solvent: {}\n'.format(T1, rgt1, slvt1)
        if rank != 1:
            data['html'] += '<br><br><u>Predicted major product (<i>p = {}</i>)</u>'.format(major_prob)
            data['html'] += '\n<br>{}'.format(major_prod)
            if major_prod != 'none found':
                url = reverse('draw_smiles', kwargs={'smiles':major_prod})
                data['html'] += '<br><img src="' + url + '">'
        elif rank == 1:
            data['html'] += '\n<br><i>Nearest neighbor got {}% yield</i>'.format(y1)

    plausible = plausible / 100.
    B = 150.
    R = 255. - (plausible > 0.5) * (plausible - 0.5) * (255. - B) * 2.
    G = 255. - (plausible < 0.5) * (0.5 - plausible) * (255. - B) * 2.
    data['html_color'] = str('#%02x%02x%02x' % (int(R), int(G), int(B)))
    
    # Save response
    DONE_SYNTH_PREDICTIONS['{}-{}'.format(smiles, synth_mincount)] = data
    return JsonResponse(data)
