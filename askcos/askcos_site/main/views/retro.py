from django.shortcuts import render, HttpResponse, redirect
from django.template.loader import render_to_string
from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.conf import settings
import django.contrib.auth.views
from pymongo.message import bson
from bson.objectid import ObjectId
from collections import defaultdict
import time
import numpy as np
import json
import os

from ..globals import RetroTransformer, RETRO_FOOTNOTE, \
    RETRO_CHIRAL_FOOTNOTE

from ..utils import ajax_error_wrapper, resolve_smiles
from .price import price_smiles_func
from .users import can_control_robot, can_avoid_banned_chemicals
from ..forms import SmilesInputForm
from ..models import BlacklistedReactions, BlacklistedChemicals

from makeit import global_config as gc
from rdkit import Chem

from askcos_site.askcos_celery.treebuilder.tb_c_worker import get_top_precursors as get_top_precursors_c
from askcos_site.askcos_celery.treebuilder.tb_worker import get_top_precursors
from askcos_site.askcos_celery.treebuilder.tb_coordinator import get_buyable_paths
from askcos_site.askcos_celery.treebuilder.tb_coordinator_mcts import get_buyable_paths as get_buyable_paths_mcts

with open(gc.BAN_LIST_PATH) as f:
    ban_list = json.load(f)
BANNED_SMILES = [smi for sublist in ban_list.values() for smi in sublist if smi is not None]
BANNED_SMILES = [
    Chem.MolToSmiles(
        Chem.MolFromSmiles(smi), True
    ) 
    for smi in BANNED_SMILES
]

def is_banned(request, smiles):
    if can_avoid_banned_chemicals(request):
        return False
    if smiles in BANNED_SMILES:
        return True
    return False

#@login_required
def retro(request, smiles=None, chiral=True, mincount=0, max_n=200):
    '''
    Retrosynthesis homepage
    '''
    context = {}

    # Set default inputs
    context['form'] = {}
    context['form']['template_prioritization'] = request.session.get('template_prioritization', 'Relevance')
    context['form']['precursor_prioritization'] = request.session.get('precursor_prioritization', 'RelevanceHeuristic')
    context['form']['template_count'] = request.session.get('template_count', '100')
    context['form']['max_cum_prob'] = request.session.get('max_cum_prob', '0.995')
    context['form']['filter_threshold'] = request.session.get('filter_threshold', '0.75')

    print(request)
    if request.method == 'POST':
        print(request)
        smiles = str(request.POST['smiles'])  # not great...
        context['form']['template_prioritization'] = str(
            request.POST['template_prioritization'])
        request.session['template_prioritization'] = context['form']['template_prioritization']
        context['form']['precursor_prioritization'] = str(
            request.POST['precursor_prioritization'])
        request.session['precursor_prioritization'] = context['form']['precursor_prioritization']

        if 'template_count' in request.POST:
            context['form']['template_count'] = str(request.POST['template_count'])
            request.session['template_count'] = context['form']['template_count']
        if 'max_cum_prob' in request.POST:
            context['form']['max_cum_prob'] = str(request.POST['max_cum_prob'])
            request.session['max_cum_prob'] = context['form']['max_cum_prob']
        if 'filter_threshold' in request.POST:
            context['form']['filter_threshold'] = str(request.POST['filter_threshold'])
            request.session['filter_threshold'] = context['form']['filter_threshold']

        smiles = resolve_smiles(smiles)
        if smiles is None:
            context['err'] = 'Could not parse!'
            return render(request, 'retro.html', context)

    if smiles is not None and not is_banned(request, smiles):

        # OLD: ALWAYS CHIRAL NOW
        # if 'retro_lit' in request.POST: return redirect('retro_lit_target', smiles=smiles)
        # if 'retro' in request.POST:
        #     return retro_target(request, smiles, chiral=False)
        # if 'retro_chiral' in request.POST:
        #     return retro_target(request, smiles, chiral=True)

        # Look up target
        smiles_img = reverse('draw_smiles', kwargs={'smiles': smiles})
        context['target'] = {
            'smiles': smiles,
            'img': smiles_img
        }

        # Perform retrosynthesis
        context['form']['smiles'] = smiles
        template_prioritization = context['form']['template_prioritization']
        precursor_prioritization = context['form']['precursor_prioritization']
        filter_threshold = context['form']['filter_threshold']

        try:
            template_count = int(context['form']['template_count'])
            if (template_count < 1):
                raise ValueError
        except ValueError:
            context['err'] = 'Invalid template count specified!'
            return render(request, 'retro.html', context)

        try:
            max_cum_prob = float(context['form']['max_cum_prob'])
            if (max_cum_prob <= 0):
                raise ValueError
        except ValueError:
            context['err'] = 'Invalid maximum cumulative probability specified!'
            return render(request, 'retro.html', context)

        try:
            filter_threshold = float(context['form']['filter_threshold'])
            filter_threshold = min(1, filter_threshold)
            filter_threshold = max(0, filter_threshold)
            apply_fast_filter = filter_threshold > 0
        except ValueError:
            context['err'] = 'Invalid filter threshold specified!'
            return render(request, 'retro.html', context)


        if template_prioritization == 'Popularity':
            template_count = 1e9

        print('Retro expansion conditions:')
        print(smiles)
        print(template_prioritization)
        print(precursor_prioritization)
        print(template_count)
        print(max_cum_prob)
        print(filter_threshold)

        startTime = time.time()
        if chiral:

            res = get_top_precursors_c.delay(
                smiles, template_prioritization, precursor_prioritization, mincount=0, max_branching=max_n,
                template_count=template_count, max_cum_prob=max_cum_prob, apply_fast_filter=apply_fast_filter, filter_threshold=filter_threshold)
            (smiles, precursors) = res.get(300)
            # allow up to 5 minutes...can be pretty slow
            context['precursors'] = precursors
            context['footnote'] = RETRO_CHIRAL_FOOTNOTE
        else:

            # Use apply_async so we can force high priority
            res = get_top_precursors.delay(smiles, template_prioritization, precursor_prioritization,
                mincount=0, max_branching=max_n, template_count=template_count, max_cum_prob=max_cum_prob, apply_fast_filter=apply_fast_filter, filter_threshold=filter_threshold)
            (smiles, precursors) = res.get(300)
            context['precursors'] = precursors
            context['footnote'] = RETRO_FOOTNOTE
        context['time'] = '%0.3f' % (time.time() - startTime)

        # Change 'tform' field to be reaction SMARTS, not ObjectID from Mongo
        # Also add up total number of examples
        for (i, precursor) in enumerate(context['precursors']):
            context['precursors'][i]['tforms'] = \
                [dict(RetroTransformer.lookup_id(ObjectId(_id)), **
                      {'id': str(_id)}) for _id in precursor['tforms']]
            context['precursors'][i]['mols'] = []
            # Overwrite num examples
            context['precursors'][i]['num_examples'] = sum(
                tform['count'] for tform in precursor['tforms'])
            for smiles in precursor['smiles_split']:
                ppg = price_smiles_func(smiles)
                context['precursors'][i]['mols'].append({
                    'smiles': smiles,
                    'ppg': '${}/g'.format(ppg) if ppg else 'cannot buy',
                })
                
    elif smiles is not None:
        context['err'] = 'ASKCOS does not provide results for compounds on restricted lists such as the CWC and DEA schedules'
    else:


        # Define suggestions
        context['suggestions'] = [
            {'name': 'Diphenhydramine', 'smiles': 'CN(C)CCOC(c1ccccc1)c2ccccc2'},
            {'name': 'Fluconazole', 'smiles': 'OC(Cn1cncn1)(Cn2cncn2)c3ccc(F)cc3F'},
            {'name': 'Nevirapine', 'smiles': 'Cc1ccnc2N(C3CC3)c4ncccc4C(=O)Nc12'},
            {'name': 'Atropine', 'smiles': 'CN1C2CCC1CC(C2)OC(=O)C(CO)c3ccccc3'},
            {'name': 'Diazepam', 'smiles': 'CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13'},
        ]


 
    return render(request, 'retro.html', context)

#@login_required
def retro_target(request, smiles):
    '''
    Given a target molecule, render page
    '''
    return retro(request, smiles=smiles)


def retro_interactive(request, target=None):
    '''Builds an interactive retrosynthesis page'''

    context = {}
    context['warn'] = 'If requests seem to take a long time, check the <a href="/status/">Server Status</a> page to see which resources are currently being used!'

    context['max_depth_default'] = 4
    context['max_branching_default'] = 20
    context['retro_mincount_default'] = settings.RETRO_TRANSFORMS['mincount']
    context['synth_mincount_default'] = settings.SYNTH_TRANSFORMS['mincount']
    context['expansion_time_default'] = 60
    context['max_ppg_default'] = 100
    context['template_count_default'] = 100
    context['template_prioritization'] = 'Relevance'
    context['max_cum_prob_default'] = 0.995
    context['precursor_prioritization'] = 'RelevanceHeuristic'
    context['forward_scorer'] = 'Template_Free'
    context['filter_threshold_default'] = 0.75

    if target is not None:
        context['target_mol'] = target

    return render(request, 'retro_interactive.html', context)


def retro_interactive_mcts(request, target=None):
    '''Builds an interactive retrosynthesis page'''

    context = {}
    context['warn'] = '<div style="text-align: center">If requests seem to take a long time, check the <a href="/status/" target="_blank">Server Status</a> page to see which resources are currently being used!</div>'

    context['max_depth_default'] = 4
    context['max_branching_default'] = 20
    context['retro_mincount_default'] = settings.RETRO_TRANSFORMS['mincount']
    context['synth_mincount_default'] = settings.SYNTH_TRANSFORMS['mincount']
    context['expansion_time_default'] = 60
    context['max_ppg_default'] = 100
    context['template_count_default'] = 100
    context['template_prioritization'] = 'Relevance'
    context['max_cum_prob_default'] = 0.995
    context['precursor_prioritization'] = 'RelevanceHeuristic'
    context['forward_scorer'] = 'Template_Free'
    context['filter_threshold_default'] = 0.75

    if target is not None:
        context['target_mol'] = target

    if request.user.is_authenticated():
        context['logged_in'] = True
    else:
        context['logged_in'] = False

    return render(request, 'retro_interactive_mcts.html', context)


@ajax_error_wrapper
def ajax_start_retro_celery(request):
    '''Start builder'''
    data = {'err': False}

    smiles = request.GET.get('smiles', None)
    max_depth = int(request.GET.get('max_depth', 4))
    max_branching = int(request.GET.get('max_branching', 25))
    retro_mincount = int(request.GET.get('retro_mincount', 0))
    expansion_time = int(request.GET.get('expansion_time', 60))
    max_ppg = int(request.GET.get('max_ppg', 10))
    chiral = json.loads(request.GET.get('chiral', 'true'))
    precursor_prioritization = request.GET.get(
        'precursor_prioritization', 'RelevanceHeuristic')
    template_prioritization = request.GET.get(
        'template_prioritization', 'Relevance')
    template_count = int(request.GET.get('template_count', '100'))
    max_cum_prob = float(request.GET.get('max_cum_prob', '0.995'))
    chemical_property_logic = str(request.GET.get('chemical_property_logic', 'none'))
    max_chemprop_c = int(request.GET.get('max_chemprop_c', '0'))
    max_chemprop_n = int(request.GET.get('max_chemprop_n', '0'))
    max_chemprop_o = int(request.GET.get('max_chemprop_o', '0'))
    max_chemprop_h = int(request.GET.get('max_chemprop_h', '0'))
    chemical_popularity_logic = str(request.GET.get('chemical_popularity_logic', 'none'))
    min_chempop_reactants = int(request.GET.get('min_chempop_reactants', 5))
    min_chempop_products = int(request.GET.get('min_chempop_products', 5))

    filter_threshold = float(request.GET.get('filter_threshold', 0.75))
    apply_fast_filter = filter_threshold > 0

    blacklisted_reactions = list(set(
        [x.smiles for x in BlacklistedReactions.objects.filter(user=request.user, active=True)]))
    forbidden_molecules = list(set(
        [x.smiles for x in BlacklistedChemicals.objects.filter(user=request.user, active=True)]))

    if template_prioritization == 'Popularity':
            template_count = 1e9

    default_val = 1e9 if chemical_property_logic == 'and' else 0
    max_natom_dict = defaultdict(lambda: default_val, {
        'logic': chemical_property_logic,
        'C': max_chemprop_c,
        'N': max_chemprop_n,
        'O': max_chemprop_o,
        'H': max_chemprop_h,
    })
    min_chemical_history_dict = {
        'logic': chemical_popularity_logic,
        'as_reactant': min_chempop_reactants,
        'as_product': min_chempop_products,
    }
    print('Tree building {} for user {} ({} forbidden reactions)'.format(
        smiles, request.user, len(blacklisted_reactions)))
    print('Using chemical property logic: {}'.format(max_natom_dict))
    print('Using chemical popularity logic: {}'.format(min_chemical_history_dict))

    res = get_buyable_paths.delay(smiles, template_prioritization, precursor_prioritization,
                                  mincount=retro_mincount, max_branching=max_branching, max_depth=max_depth,
                                  max_ppg=max_ppg, max_time=expansion_time, max_trees=500, reporting_freq=5,
                                  chiral=chiral, known_bad_reactions=blacklisted_reactions,
                                  forbidden_molecules=forbidden_molecules,
                                  max_cum_template_prob=max_cum_prob, template_count=template_count,
                                  max_natom_dict=max_natom_dict, min_chemical_history_dict=min_chemical_history_dict,
                                  apply_fast_filter=apply_fast_filter, filter_threshold=filter_threshold)
    (tree_status, trees) = res.get(expansion_time * 3)

    # print(trees)

    (num_chemicals, num_reactions, at_depth) = tree_status
    data['html_stats'] = 'After expanding (with {} banned reactions, {} banned chemicals), {} total chemicals and {} total reactions'.format(
        len(blacklisted_reactions), len(forbidden_molecules), num_chemicals, num_reactions)
    for (depth, count) in sorted(at_depth.items(), key=lambda x: x[0]):
        label = 'Could not format label...?'
        if int(float(depth)) == float(depth):
            label = 'chemicals'
        else:
            label = 'reactions'
        data[
            'html_stats'] += '<br>   at depth {}, {} {}'.format(depth, count, label)

    if trees:
        data['html_trees'] = render_to_string('trees_only.html',
                                              {'trees': trees, 'can_control_robot': can_control_robot(request)})
    else:
        data['html_trees'] = render_to_string('trees_none.html', {})

    # Save to session in case user wants to export
    request.session['last_retro_interactive'] = trees
    print('Saved {} trees to {} session'.format(
        len(trees), request.user.get_username()))

    return JsonResponse(data)


@ajax_error_wrapper
def ajax_start_retro_mcts_celery(request):
    '''Start builder'''
    data = {'err': False}

    smiles = request.GET.get('smiles', None)
    
    if is_banned(request, smiles):
        data['html_trees'] = 'ASKCOS does not provide results for compounds on restricted lists such as the CWC and DEA schedules'
        return JsonResponse(data)
    
    max_depth = int(request.GET.get('max_depth', 4))
    max_branching = int(request.GET.get('max_branching', 25))
    expansion_time = int(request.GET.get('expansion_time', 60))
    max_ppg = int(request.GET.get('max_ppg', 10))
    template_count = int(request.GET.get('template_count', '100'))
    max_cum_prob = float(request.GET.get('max_cum_prob', '0.995'))
    chemical_property_logic = str(request.GET.get('chemical_property_logic', 'none'))
    max_chemprop_c = int(request.GET.get('max_chemprop_c', '0'))
    max_chemprop_n = int(request.GET.get('max_chemprop_n', '0'))
    max_chemprop_o = int(request.GET.get('max_chemprop_o', '0'))
    max_chemprop_h = int(request.GET.get('max_chemprop_h', '0'))
    chemical_popularity_logic = str(request.GET.get('chemical_popularity_logic', 'none'))
    min_chempop_reactants = int(request.GET.get('min_chempop_reactants', 5))
    min_chempop_products = int(request.GET.get('min_chempop_products', 5))
    filter_threshold = float(request.GET.get('filter_threshold', 0.75))
    apply_fast_filter = filter_threshold > 0
    return_first = json.loads(request.GET.get('return_first', 'false'))

    if request.user.is_authenticated():
        blacklisted_reactions = list(set(
            [x.smiles for x in BlacklistedReactions.objects.filter(user=request.user, active=True)]))
        forbidden_molecules = list(set(
            [x.smiles for x in BlacklistedChemicals.objects.filter(user=request.user, active=True)]))
    else:
        blacklisted_reactions = []
        forbidden_molecules = []

    default_val = 1e9 if chemical_property_logic == 'and' else 0
    max_natom_dict = defaultdict(lambda: default_val, {
        'logic': chemical_property_logic,
        'C': max_chemprop_c,
        'N': max_chemprop_n,
        'O': max_chemprop_o,
        'H': max_chemprop_h,
    })
    min_chemical_history_dict = {
        'logic': chemical_popularity_logic,
        'as_reactant': min_chempop_reactants,
        'as_product': min_chempop_products,
    }
    print('Tree building {} for user {} ({} forbidden reactions)'.format(
        smiles, request.user.id, len(blacklisted_reactions)))
    print('Using chemical property logic: {}'.format(max_natom_dict))
    print('Using chemical popularity logic: {}'.format(min_chemical_history_dict))
    print('Returning as soon as any pathway found? {}'.format(return_first))

    res = get_buyable_paths_mcts.delay(smiles, max_branching=max_branching, max_depth=max_depth,
                                  max_ppg=max_ppg, expansion_time=expansion_time, max_trees=500,
                                  known_bad_reactions=blacklisted_reactions,
                                  forbidden_molecules=forbidden_molecules,
                                  max_cum_template_prob=max_cum_prob, template_count=template_count,
                                  max_natom_dict=max_natom_dict, min_chemical_history_dict=min_chemical_history_dict,
                                  apply_fast_filter=apply_fast_filter, filter_threshold=filter_threshold,
                                  return_first=return_first)
    (tree_status, trees) = res.get(expansion_time * 3)


    # print(trees)

    (num_chemicals, num_reactions, _) = tree_status
    data['html_stats'] = 'After expanding (with {} banned reactions, {} banned chemicals), {} total chemicals and {} total reactions'.format(
        len(blacklisted_reactions), len(forbidden_molecules), num_chemicals, num_reactions)

    if trees:
        data['html_trees'] = render_to_string('trees_only.html',
                                              {'trees': trees, 'can_control_robot': can_control_robot(request)})
    else:
        data['html_trees'] = render_to_string('trees_none.html', {})

    # Save to session in case user wants to export
    request.session['last_retro_interactive'] = trees
    # print('Saved {} trees to {} session'.format(
    #     len(trees), request.user.get_username()))

    return JsonResponse(data)
