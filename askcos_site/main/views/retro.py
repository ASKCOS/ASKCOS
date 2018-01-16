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

from ..globals import RetroTransformer, RETRO_FOOTNOTE, \
    RETRO_CHIRAL_FOOTNOTE

from ..utils import ajax_error_wrapper, resolve_smiles
from .price import price_smiles_func
from .users import can_control_robot
from ..forms import SmilesInputForm
from ..models import BlacklistedReactions

@login_required
def retro(request, smiles=None, chiral=True, mincount=0, max_n=200):
    '''
    Retrosynthesis homepage
    '''
    context = {}

    # Set default inputs
    context['form'] = {}
    context['form']['template_prioritization'] = 'Popularity'
    context['form']['precursor_prioritization'] = 'Heuristic'
    print(request.method)
    if request.method == 'POST':
        print(context['form'])
        smiles = str(request.POST['smiles']) # not great...
        context['form']['template_prioritization'] = str(request.POST['template_prioritization'])
        context['form']['precursor_prioritization'] = str(request.POST['precursor_prioritization'])

        smiles = resolve_smiles(smiles)
        if smiles is None:
            context['err'] = 'Could not parse!'
            return render(request, 'retro.html', context)

    if smiles is not None:

        # OLD: ALWAYS CHIRAL NOW
        # if 'retro_lit' in request.POST: return redirect('retro_lit_target', smiles=smiles)
        # if 'retro' in request.POST:
        #     return retro_target(request, smiles, chiral=False)
        # if 'retro_chiral' in request.POST:
        #     return retro_target(request, smiles, chiral=True)
        
        # Look up target
        smiles_img = reverse('draw_smiles', kwargs={'smiles':smiles})
        context['target'] = {
            'smiles': smiles,
            'img': smiles_img
        }

        # Perform retrosynthesis
        context['form']['smiles'] = smiles
        template_prioritization = context['form']['template_prioritization']
        precursor_prioritization = context['form']['precursor_prioritization']
        
        print('Retro expansion conditions:')
        print(smiles)
        print(template_prioritization)
        print(precursor_prioritization)

        startTime = time.time()
        if chiral:
            from askcos_site.askcos_celery.treebuilder.tb_c_worker import get_top_precursors
            res = get_top_precursors.delay(smiles, template_prioritization, precursor_prioritization, mincount=0, max_branching=max_n)
            (smiles, precursors) = res.get(300)
            context['precursors'] = precursors # allow up to 5 minutes...can be pretty slow
            context['footnote'] = RETRO_CHIRAL_FOOTNOTE
        else:
            from askcos_site.askcos_celery.treebuilder.tb_worker import get_top_precursors
            # Use apply_async so we can force high priority 
            res = get_top_precursors.apply_async(args=(smiles, template_prioritization, precursor_prioritization), 
                kwargs={'mincount':0, 'max_branching':max_n, 'raw_results':True}, 
                priority=255)
            context['precursors'] = res.get(120)
            context['footnote'] = RETRO_FOOTNOTE
        context['time'] = '%0.3f' % (time.time() - startTime)

        # Change 'tform' field to be reaction SMARTS, not ObjectID from Mongo
        # Also add up total number of examples
        for (i, precursor) in enumerate(context['precursors']):
            context['precursors'][i]['tforms'] = \
                [dict(RetroTransformer.lookup_id(ObjectId(_id)), **{'id':str(_id)}) for _id in precursor['tforms']]
            context['precursors'][i]['mols'] = []
            # Overwrite num examples
            context['precursors'][i]['num_examples'] = sum(tform['count'] for tform in precursor['tforms'])
            for smiles in precursor['smiles_split']:
                ppg = price_smiles_func(smiles)
                context['precursors'][i]['mols'].append({
                    'smiles': smiles,
                    'ppg': '${}/g'.format(ppg) if ppg else 'cannot buy'
                })

    else:

        # Define suggestions
        context['suggestions'] = [
            {'name': 'Diphenhydramine', 'smiles': 'CN(C)CCOC(c1ccccc1)c2ccccc2'},
            {'name': 'Fluconazole', 'smiles': 'OC(Cn1cncn1)(Cn2cncn2)c3ccc(F)cc3F'},
            {'name': 'Nevirapine', 'smiles': 'Cc1ccnc2N(C3CC3)c4ncccc4C(=O)Nc12'},
            {'name': 'Atropine', 'smiles': 'CN1C2CCC1CC(C2)OC(=O)C(CO)c3ccccc3'},
            {'name': 'Diazepam', 'smiles': 'CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13'},
            {'name': 'Hydroxychloroquine', 'smiles': 'CCN(CCO)CCCC(C)Nc1ccnc2cc(Cl)ccc12'},
            {'name': 'Ibuprofen', 'smiles': 'CC(C)Cc1ccc(cc1)C(C)C(O)=O'},
            {'name': 'Tramadol', 'smiles': 'CN(C)C[C@H]1CCCC[C@@]1(C2=CC(=CC=C2)OC)O'},
            {'name': 'Lamivudine', 'smiles': 'NC1=NC(=O)N(C=C1)[C@@H]2CS[C@H](CO)O2'},
            {'name': 'Pregabalin', 'smiles': 'CC(C)C[C@H](CN)CC(O)=O'},
            {'name': 'Naproxen', 'smiles': 'COc1ccc2cc(ccc2c1)C(C)C(O)=O'},
            {'name': 'Imatinib', 'smiles': 'CN1CCN(CC1)Cc2ccc(cc2)C(=O)Nc3ccc(C)c(Nc4nccc(n4)c5cccnc5)c3'},
            {'name': 'Quinapril', 'smiles': 'CCOC(=O)[C@H](CCc1ccccc1)N[C@@H](C)C(=O)N2Cc3ccccc3C[C@H]2C(O)=O'},
            {'name': 'Atorvastatin', 'smiles': 'CC(C)c1n(CC[C@@H](O)C[C@@H](O)CC(O)=O)c(c2ccc(F)cc2)c(c3ccccc3)c1C(=O)Nc4ccccc4'},
            {'name': 'Bortezomib', 'smiles': 'CC(C)C[C@@H](NC(=O)[C@@H](Cc1ccccc1)NC(=O)c2cnccn2)B(O)O'},
            {'name': 'Itraconazole', 'smiles': 'CCC(C)N1N=CN(C1=O)c2ccc(cc2)N3CCN(CC3)c4ccc(OC[C@H]5CO[C@@](Cn6cncn6)(O5)c7ccc(Cl)cc7Cl)cc4'},
            {'name': '6-Carboxytetramethylrhodamine', 'smiles': 'CN(C)C1=CC2=C(C=C1)C(=C3C=CC(=[N+](C)C)C=C3O2)C4=C(C=CC(=C4)C(=O)[O-])C(=O)O'},
            {'name': '6-Carboxytetramethylrhodamine isomer', 'smiles': 'CN(C)c1ccc2c(c1)Oc1cc(N(C)C)ccc1C21OC(=O)c2ccc(C(=O)O)c1c2'},
            {'name': '(S)-Warfarin', 'smiles': 'CC(=O)C[C@@H](C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O'},
            {'name': 'Tranexamic Acid', 'smiles': 'NC[C@@H]1CC[C@H](CC1)C(O)=O'},
        ]

    context['footnote'] = RETRO_CHIRAL_FOOTNOTE
    return render(request, 'retro.html', context)

@login_required
def retro_target(request, smiles):
    '''
    Given a target molecule, render page
    '''
    return retro(request, smiles=smiles)

    # Render form with target
    context = {}
    smiles = resolve_smiles(smiles)
    context['form'] = SmilesInputForm({'smiles': smiles})
    if smiles is None:
        context['err'] = 'Could not parse!'
        return render(request, 'retro.html', context)

    # Look up target
    smiles_img = reverse('draw_smiles', kwargs={'smiles':smiles})
    context['target'] = {
        'smiles': smiles,
        'img': smiles_img
    }

    # Perform retrosynthesis
    template_prioritization = 'Popularity'
    precursor_prioritization = 'Heuristic'
    startTime = time.time()
    if chiral:
        from askcos_site.askcos_celery.treebuilder.tb_c_worker import get_top_precursors
        res = get_top_precursors.delay(smiles, template_prioritization, precursor_prioritization, mincount=0, max_branching=max_n)
        (smiles, precursors) = res.get(300)
        context['precursors'] = precursors # allow up to 5 minutes...can be pretty slow
        context['footnote'] = RETRO_CHIRAL_FOOTNOTE
    else:
        from askcos_site.askcos_celery.treebuilder.tb_worker import get_top_precursors
        # Use apply_async so we can force high priority 
        res = get_top_precursors.apply_async(args=(smiles, template_prioritization, precursor_prioritization), 
            kwargs={'mincount':0, 'max_branching':max_n, 'raw_results':True}, 
            priority=255)
        context['precursors'] = res.get(120)
        context['footnote'] = RETRO_FOOTNOTE
    context['time'] = '%0.3f' % (time.time() - startTime)

    # Change 'tform' field to be reaction SMARTS, not ObjectID from Mongo
    # Also add up total number of examples
    for (i, precursor) in enumerate(context['precursors']):
        context['precursors'][i]['tforms'] = \
            [dict(RetroTransformer.lookup_id(ObjectId(_id)), **{'id':str(_id)}) for _id in precursor['tforms']]
        context['precursors'][i]['mols'] = []
        # Overwrite num examples
        context['precursors'][i]['num_examples'] = sum(tform['count'] for tform in precursor['tforms'])
        for smiles in precursor['smiles_split']:
            ppg = price_smiles_func(smiles)
            context['precursors'][i]['mols'].append({
                'smiles': smiles,
                'ppg': '${}/g'.format(ppg) if ppg else 'cannot buy'
            })

    return render(request, 'retro.html', context)

@login_required
def retro_interactive(request, target=None):
    '''Builds an interactive retrosynthesis page'''

    context = {}
    context['warn'] = 'The worker pool is not set up for autoscaling; there is a chance that all of the tree building coordinators and workers will be occupied when you try to run a target.'

    context['max_depth_default'] = 4
    context['max_branching_default'] = 20
    context['retro_mincount_default'] = settings.RETRO_TRANSFORMS['mincount']
    context['synth_mincount_default'] = settings.SYNTH_TRANSFORMS['mincount']
    context['expansion_time_default'] = 60
    context['max_ppg_default'] = 100

    if target is not None:
        context['target_mol'] = target

    return render(request, 'retro_interactive.html', context)


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
    chiral = json.loads(request.GET.get('chiral', 'false'))

    blacklisted_reactions = list(set([x.smiles for x in BlacklistedReactions.objects.filter(user=request.user, active=True)]))

    from askcos_site.askcos_celery.treebuilder.coordinator import get_buyable_paths
    res = get_buyable_paths.delay(smiles, mincount=retro_mincount, max_branching=max_branching, max_depth=max_depth, 
        max_ppg=max_ppg, max_time=expansion_time, max_trees=500, reporting_freq=5, chiral=chiral,
        known_bad_reactions=blacklisted_reactions)
    (tree_status, trees) = res.get(expansion_time * 3)
    print('Tree building {} for user {} ({} forbidden reactions)'.format(
        smiles, request.user, len(blacklisted_reactions)))
    # print(trees)

    (num_chemicals, num_reactions, at_depth) = tree_status
    data['html_stats'] = 'After expanding, {} total chemicals and {} total reactions'.format(num_chemicals, num_reactions)
    for (depth, count) in sorted(at_depth.iteritems(), key=lambda x: x[0]):
        label = 'Could not format label...?'
        if int(float(depth)) == float(depth):
            label = 'chemicals'
        else:
            label = 'reactions'
        data['html_stats'] += '<br>   at depth {}, {} {}'.format(depth, count, label)

    if trees:
        data['html_trees'] = render_to_string('trees_only.html', 
            {'trees': trees, 'can_control_robot': can_control_robot(request)})
    else:
        data['html_trees'] = 'No trees resulting in buyable chemicals found! If the program is having trouble with your target, you may want to explore the One-Step Retrosynthesis options and help guide the search.'
    
    # Save to session in case user wants to export
    request.session['last_retro_interactive'] = trees
    print('Saved {} trees to {} session'.format(len(trees), request.user.get_username()))

    return JsonResponse(data)     


