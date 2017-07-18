from django.shortcuts import render, HttpResponse, redirect
from django.template.loader import render_to_string
from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.conf import settings
import django.contrib.auth.views
from forms import SmilesInputForm, DrawingInputForm, is_valid_smiles
from pymongo.message import bson
from bson.objectid import ObjectId
from datetime import datetime
import time
import numpy as np 
import json
import os

import rdkit.Chem as Chem 
import urllib2

from askcos_site.askcos_celery.chiralretro.coordinator import get_top_chiral_precursors
from askcos_site.askcos_celery.treebuilder.worker import get_top_precursors

from askcos_site.main.globals import RetroTransformer, RETRO_FOOTNOTE, \
    SynthTransformer, SYNTH_FOOTNOTE, REACTION_DB, INSTANCE_DB, CHEMICAL_DB, \
    BUYABLE_DB, SOLVENT_DB, Pricer, TransformerOnlyKnown, \
    PREDICTOR_FOOTNOTE, DONE_SYNTH_PREDICTIONS, RETRO_CHIRAL_FOOTNOTE


def log_this_request(method):
    def f(*args, **kwargs):
        try:
            print('User %s requested view %s with args %r and kwargs %r' % \
                args[0].user.get_username(), method.__name__, args, kwargs)
        except Exception as e:
            print(e)
        return method(*args, **kwargs)
    return f

@login_required
def index(request):
    '''
    Homepage
    '''
    print('Somebody loaded the index page!')
    return render(request, 'index.html')

def login(request):
    '''
    User login
    '''
    return django.contrib.auth.views.login(request, template_name = 'login.html')

def logout(request):
    '''
    User logout
    '''
    return django.contrib.auth.views.logout(request, template_name = 'logout.html')

from models import SavedResults
@login_required
def user_saved_results(request, err=None):
    saved_results = SavedResults.objects.filter(user=request.user)
    return render(request, 'saved_results.html', {'saved_results':saved_results, 'err': err})

@login_required
def user_saved_results_id(request, _id=-1):
    saved_result = SavedResults.objects.filter(user=request.user, id=_id)
    if saved_result.count() == 0:
        return user_saved_results(request, err='Could not find that ID')
    with open(saved_result[0].fpath, 'r') as fid:
        html = fid.read()
    return render(request, 'saved_results_id.html', 
        {'saved_result':saved_result[0], 'html':html})

@login_required
def user_saved_results_del(request, _id=-1):
    SavedResults.objects.filter(user=request.user, id=_id).delete()
    return user_saved_results(request, err=None)

@login_required
def ajax_user_save_page(request):
    html = request.GET.get('html', None)
    if html is None:
        data = {'err': 'Could not get HTML to save'}
        return JsonResponse(data)
    print('Got request to save a page')
    now = datetime.now()
    unique_str = '%i.txt' % hash((now, request.user))
    fpath = os.path.join(settings.LOCAL_STORAGE['user_saves'], unique_str)
    _id = SavedResults.objects.create(user=request.user, 
        description='',
        created=datetime.now(),
        fpath=fpath)
    print('Created saved object {}'.format(_id))
    with open(fpath, 'w') as fid:
        fid.write(html)
    print('Wrote to {}'.format(fpath))

    return JsonResponse({'err': False})

@login_required
def retro(request):
    '''
    Retrosynthesis homepage
    '''
    context = {}

    if request.method == 'POST':
        context['form'] = SmilesInputForm(request.POST)
        if not context['form'].is_valid():
            context['err'] = 'Could not parse!'
        else:
            # Identify target
            smiles = context['form'].cleaned_data['smiles']
            if 'retro_lit' in request.POST: return redirect('retro_lit_target', smiles=smiles)
            if 'retro' in request.POST:
                return retro_target(request, smiles, chiral=False)
            if 'retro_chiral' in request.POST:
                return retro_target(request, smiles, chiral=True)
    else:
        context['form'] = SmilesInputForm()

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
        {'name': '(S)-Warfarin', 'smiles': 'CC(=O)C[C@@H](C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O'},
        {'name': 'Tranexamic Acid', 'smiles': 'NC[C@@H]1CC[C@H](CC1)C(O)=O'},
    ]

    context['footnote'] = RETRO_CHIRAL_FOOTNOTE
    return render(request, 'retro.html', context)

@login_required
def retro_target(request, smiles, chiral=True, max_n=200):
    '''
    Given a target molecule, render page
    '''

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
    startTime = time.time()
    if chiral:
        res = get_top_chiral_precursors.delay(smiles, mincount=0, max_branching=max_n, raw_results=True)
        context['precursors'] = res.get(300) # allow up to 5 minutes...can be pretty slow
        context['footnote'] = RETRO_CHIRAL_FOOTNOTE
    else:
        # Use apply_async so we can force high priority 
        res = get_top_precursors.apply_async(args=(smiles,), 
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
            ppg = Pricer.lookup_smiles(smiles, alreadyCanonical=True)
            context['precursors'][i]['mols'].append({
                'smiles': smiles,
                'ppg': '${}/g'.format(ppg) if ppg else 'cannot buy'
            })

    return render(request, 'retro.html', context)

@login_required
def retro_lit_target(request, smiles, max_n = 50):
    '''
    Given a target molecule, render page
    '''

    # Render form with target
    context = {}

    # Perform retrosynthesis
    result = TransformerOnlyKnown.perform_retro(smiles)
    if result:
        context['precursors'] = result.return_top(n = 50)
        # Erase 'tform' field - we have specific RX IDs, not templates
        # Also add up total number of examples
        for (i, precursor) in enumerate(context['precursors']):
            context['precursors'][i]['rxid'] = context['precursors'][i]['tforms'][0]
            del context['precursors'][i]['tforms']
            context['precursors'][i]['mols'] = []
            for smiles in precursor['smiles_split']:
                ppg = Pricer.lookup_smiles(smiles, alreadyCanonical = True)
                context['precursors'][i]['mols'].append({
                    'smiles': smiles,
                    'ppg': '${}/g'.format(ppg) if ppg else 'cannot buy'
                })
        smiles = result.target_smiles
    else:
        context['err'] = 'Could not find a database match for {}'.format(smiles)

    context['form'] = SmilesInputForm({'smiles': smiles})
    # Look up target
    smiles_img = reverse('draw_smiles', kwargs={'smiles':smiles})
    context['target'] = {
        'smiles': smiles,
        'img': smiles_img
    }
    context['lit_only'] = True
    context['footnote'] = RETRO_LIT_FOOTNOTE
    return render(request, 'retro.html', context)

@login_required
def retro_interactive(request):
    '''Builds an interactive retrosynthesis page'''

    context = {}
    context['warn'] = 'The worker pool is not set up for autoscaling; there is a chance that all of the tree building coordinators and workers will be occupied when you try to run a target.'

    context['max_depth_default'] = 4
    context['max_branching_default'] = 20
    context['retro_mincount_default'] = settings.RETRO_TRANSFORMS['mincount']
    context['synth_mincount_default'] = settings.SYNTH_TRANSFORMS['mincount']
    context['expansion_time_default'] = 60
    context['max_ppg_default'] = 100

    return render(request, 'retro_interactive.html', context)

@login_required
def synth_interactive(request, reactants='', reagents='', solvent='toluene', 
        temperature='20', mincount='', product=None):
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
        from askcos_site.askcos_celery.contextrecommender.worker import get_context_recommendation
        res = get_context_recommendation.delay(smiles, n=1)
        contexts = res.get(60)
        if contexts is None:
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

def ajax_error_wrapper(ajax_func):
    if not settings.DEBUG:
        def ajax_func_call(*args, **kwargs):
            try:
                return ajax_func(*args, **kwargs)
            except Exception as e:
                return JsonResponse({'err':True, 'message': str(e)})
    else:
        def ajax_func_call(*args, **kwargs):
            try:
                return ajax_func(*args, **kwargs)
            except Exception as e:
                print(e)
                raise(e)

    return ajax_func_call

def resolve_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        # Try to resolve using NIH
        new_smiles = []
        for smiles in smiles.split(' and '):
            try:
                smiles = urllib2.urlopen('https://cactus.nci.nih.gov/chemical/structure/{}/smiles'.format(smiles)).read()
            except urllib2.HTTPError:
                return None
            mol = Chem.MolFromSmiles(smiles)
            if not mol: return None
            new_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))
        return '.'.join(new_smiles)
    return Chem.MolToSmiles(mol, isomericSmiles=True)

@ajax_error_wrapper
def ajax_smiles_to_image(request):
    '''Takes an Ajax call with a smiles string
    and returns the HTML for embedding an image'''

    smiles = request.GET.get('smiles', None)
    print('SMILES from Ajax: {}'.format(smiles))
    smiles = resolve_smiles(smiles)
    if smiles is None:
        return JsonResponse({'err': True})
    print('Resolved smiles -> {}'.format(smiles))

    url = reverse('draw_smiles', kwargs={'smiles':smiles})
    data = {
        'html': '<img src="' + url + '">',
        'smiles': smiles,
    }

    return JsonResponse(data)

@ajax_error_wrapper
def ajax_rxn_to_image(request):
    '''Takes an Ajax call with a rxn smiles string
    and returns the HTML for embedding an image'''

    reactants = request.GET.get('reactants', '')
    product = request.GET.get('product', '')

    reactants = resolve_smiles(reactants)
    product = resolve_smiles(product)
    smiles = reactants + '>>' + product
    print('RXN SMILES from Ajax: {}'.format(smiles))
    url = reverse('draw_reaction', kwargs={'smiles':smiles})
    data = {
        'html': '<img src="' + url + '">',
        'smiles': smiles,
        'reactants': reactants,
        'product': product,
    }
    return JsonResponse(data)

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
    print('Conditions for forward synthesis:')
    print('reactants: {}'.format(reactants))
    print('solvent: {}'.format(solvent))
    print('temp: {}'.format(temperature))
    print('reagents: {}'.format(reagents))
    print('mincount: {}'.format(mincount))

    startTime = time.time()
    from askcos_site.askcos_celery.forwardpredictor.coordinator import get_outcomes
    res = get_outcomes.delay(reactants, 
        contexts=[(temperature, reagents, solvent)], 
        mincount=mincount,
        top_n=25)
    outcomes = res.get(300)[0]

    print('Got top outcomes, length {}'.format(len(outcomes)))
    data['html_time'] = '{:.3f} seconds elapsed'.format(time.time() - startTime)

    if outcomes:
        data['html'] = render_to_string('synth_outcomes_only.html', {'outcomes': outcomes})
    else:
        data['html'] = 'No outcomes found? That is weird...'
    return JsonResponse(data)

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

    from askcos_site.askcos_celery.treebuilder.coordinator import get_buyable_paths
    res = get_buyable_paths.delay(smiles, mincount=retro_mincount, max_branching=max_branching, max_depth=max_depth, 
        max_ppg=max_ppg, max_time=expansion_time, max_trees=500, reporting_freq=5, chiral=chiral)
    (tree_status, trees) = res.get(expansion_time * 3)
    print(tree_status)
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
        data['html_trees'] = render_to_string('trees_only.html', {'trees': trees})
    else:
        data['html_trees'] = 'No trees resulting in buyable chemicals found!'
    return JsonResponse(data)     

@login_required
def evaluate_rxnsmiles(request):
    return render(request, 'evaluate.html', {})

# Clean up
def fix_rgt_cat_slvt(rgt1, cat1, slvt1):
    # Merge cat and reagent for forward predictor
    if rgt1 and cat1:
        rgt1 = rgt1 + '.' + cat1
    elif cat1:
        rgt1 = cat1
    # Reduce solvent to single one
    if '.' in slvt1:
        slvt1 = slvt1.split('.')[0]
    return (rgt1, cat1, slvt1)
def trim_trailing_period(txt):
    if txt:
        if txt[-1] == '.':
            return txt[:-1]
    return txt

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
    from askcos_site.askcos_celery.contextrecommender.worker import get_context_recommendation
    res = get_context_recommendation.delay(smiles, n=num_contexts)
    contexts = res.get(60)
    print('Got context(s)')
    print(contexts)
    if contexts is None:
        raise ValueError('Context recommender was unable to get valid context(?)')

    contexts_for_predictor = []
    for (T1, slvt1, rgt1, cat1, t1, y1) in contexts:
        slvt1 = trim_trailing_period(slvt1)
        rgt1 = trim_trailing_period(rgt1)
        cat1 = trim_trailing_period(cat1)
        (rgt1, cat1, slvt1) = fix_rgt_cat_slvt(rgt1, cat1, slvt1)
        contexts_for_predictor.append((T1, rgt1, slvt1))
    print('Cleaned contexts')

    # Run
    reactant_smiles = smiles.split('>>')[0]
    print('Running forward evaluator on {}'.format(reactant_smiles))
    if necessary_reagent and contexts_for_predictor[0][1]:
        reactant_smiles += contexts_for_predictor[0][1] # add rgt
    from askcos_site.askcos_celery.forwardpredictor.coordinator import get_outcomes
    res = get_outcomes.delay(reactant_smiles, contexts=contexts_for_predictor, mincount=synth_mincount, top_n=10)
    all_outcomes = res.get(300)
    if all([len(outcome) == 0 for outcome in all_outcomes]):
        if not verbose:
            data['html'] = 'Could not get outcomes - recommended context(s) unparseable'
            for i, (T, rgt, slvt) in enumerate(contexts_for_predictor):
                data['html'] += '<br>{}) T={}, rgt={}, slvt={}'.format(i+1, T, rgt, slvt)
            data['html_color'] = str('#%02x%02x%02x' % (int(255), int(0), int(0)))
            return JsonResponse(data)
        else:
            # TODO: expand
            data['html'] = '<h3>Could not get outcomes - recommended context(s) unparseable</h3>\n<ol>\n'
            for i, (T, rgt, slvt) in enumerate(contexts_for_predictor):
                data['html'] += '<li>Temp: {} C<br>Reagents: {}<br>Solvent: {}</li>\n'.format(T, rgt, slvt)
            data['html'] += '</ol>'
            data['html_color'] = str('#%02x%02x%02x' % (int(255), int(0), int(0)))
            return JsonResponse(data)
    plausible = [0. for i in range(len(all_outcomes))]
    ranks = ['>10' for i in range(len(all_outcomes))]
    major_prods = ['none found' for i in range(len(all_outcomes))]
    major_probs = ['n/a' for i in range(len(all_outcomes))]
    for i, outcomes in enumerate(all_outcomes):
        if len(outcomes) != 0:
            major_prods[i] = outcomes[0]['smiles']
            major_probs[i] = outcomes[0]['prob']
        for j, outcome in enumerate(outcomes):
            if outcome['smiles'] == products[0]:
                plausible[i] = float(outcome['prob'])
                ranks[i] = j + 1
                break
    best_context_i = np.argmax(plausible)
    plausible = plausible[best_context_i]
    rank = ranks[best_context_i]
    best_context = contexts_for_predictor[best_context_i]
    major_prod = major_prods[best_context_i]
    major_prob = major_probs[best_context_i]

    # Report
    print('Recommended context(s): {}'.format(best_context))
    print('Plausibility: {}'.format(plausible))
    (T1, rgt1, slvt1) = best_context

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


    B = 150.
    R = 255. - (plausible > 0.5) * (plausible - 0.5) * (255. - B) * 2.
    G = 255. - (plausible < 0.5) * (0.5 - plausible) * (255. - B) * 2.
    data['html_color'] = str('#%02x%02x%02x' % (int(R), int(G), int(B)))
    
    # Save response
    DONE_SYNTH_PREDICTIONS['{}-{}'.format(smiles, synth_mincount)] = data
    return JsonResponse(data)

# @login_required
# def synth(request):
#     '''
#     Forward synthesis homepage
#     '''
#     context = {}

#     if request.method == 'POST':
#         context['form'] = SmilesInputForm(request.POST)
#         if not context['form'].is_valid():
#             context['err'] = 'Could not parse!'
#         else:
#             # Identify target
#             smiles = context['form'].cleaned_data['smiles']
#             if is_valid_smiles(smiles):
#                 return redirect('synth_target', smiles = smiles)
#             else:
#                 context['err'] = 'Invalid SMILES string: {}'.format(smiles)
#     else:
#         context['form'] = SmilesInputForm()

#     context['footnote'] = SYNTH_FOOTNOTE
#     return render(request, 'synth.html', context)

# @login_required
# def synth_target(request, smiles, max_n = 50):
#     '''
#     Given a set of reactants as a single SMILES string, find products
#     '''

#     # Render form with target
#     context = {}
#     context['form'] = SmilesInputForm({'smiles': smiles})

#     # Look up target
#     smiles_img = reverse('draw_smiles', kwargs={'smiles':smiles})
#     context['target'] = {
#         'smiles': smiles,
#         'img': smiles_img,
#     }

#     # Perform forward synthesis
#     result = SynthTransformer.perform_forward(smiles)
#     context['products'] = result.return_top(n = 50)
#     #print(context['products'])
#     # Change 'tform' field to be reaction SMARTS, not ObjectID from Mongo
#     # (add "id" field to be string version of ObjectId "_id")
#     for (i, product) in enumerate(context['products']):
#         context['products'][i]['tforms'] = \
#             [dict(SynthTransformer.lookup_id(_id), **{'id':str(_id)}) for _id in product['tforms']]

#     context['footnote'] = SYNTH_FOOTNOTE
#     return render(request, 'synth.html', context)

def fancyjoin(lst, nonemessage='(none)'):
    if not lst:
        return nonemessage
    if len(lst) == 1:
        return lst[0]
    if len(lst) == 2:
        return '%s and %s' % (lst[0], lst[1])
    return ', '.join(lst[:-1]) + ', and %s' % lst[-1]

@login_required
def template_target(request, id):
    '''
    Examines a template from its id in the database
    where id == str(_id)

    Reaxys templates refer to instances, but we need the
    reactions for visualization
    '''
    context = {}

    transform = RetroTransformer.lookup_id(ObjectId(id))
    if not transform: transform = SynthTransformer.lookup_id(ObjectId(id))

    if not transform:
        context['err'] = 'Transform not found'
        return render(request, 'template.html', context)
    context['template'] = transform
    reference_ids = transform['references']

    references = []
    rx_docs = {}; xrn_to_smiles = {}
    context['total_references'] = len(reference_ids)
    for rxd_doc in INSTANCE_DB.find({'_id': {'$in': reference_ids}}):
        rx_id = int(rxd_doc['_id'].split('-')[0])
        rx_doc = REACTION_DB.find_one({'_id': rx_id})
        if rx_doc is None: rx_doc = {}
        ref = {
            'rx_id': rx_id,
            'rxd_num': rxd_doc['_id'].split('-')[1],
            'label': rxd_doc['_id'],
            'T': rxd_doc['RXD_T'] if rxd_doc['RXD_T'] != -1 else 'unk',
            'y': float(rxd_doc['RXD_NYD']) if rxd_doc['RXD_NYD'] != -1 else 'unk',
            't': rxd_doc['RXD_TIM'] if rxd_doc['RXD_TIM'] != -1 else 'unk',
            'smiles': rx_doc.get('RXN_SMILES', 'no smiles found'),
            'nvar': rx_doc.get('RX_NVAR', '?'),
            'cond': fancyjoin(rxd_doc['RXD_COND'], nonemessage=''),
            'ded': rxd_doc['RXD_DED'][0],
        }

        def rxn_lst_to_name_lst(xrn_lst):
            lst = []
            for xrn in xrn_lst:
                if xrn not in xrn_to_smiles: 
                    chem_doc = CHEMICAL_DB.find_one({'_id': xrn})
                    if 'IDE_CN' not in chem_doc:
                        if 'SMILES' not in chem_doc:
                            xrn_to_smiles[xrn] = 'Chem-%i' % xrn
                        else:
                            xrn_to_smiles[xrn] = chem_doc['SMILES']
                    else:
                        xrn_to_smiles[xrn] = chem_doc['IDE_CN']
                lst.append(xrn_to_smiles[xrn])
            return lst

        ref['reagents'] = fancyjoin(rxn_lst_to_name_lst(rxd_doc['RXD_RGTXRN']))
        ref['solvents'] = fancyjoin(rxn_lst_to_name_lst(rxd_doc['RXD_SOLXRN']))
        ref['catalysts'] = fancyjoin(rxn_lst_to_name_lst(rxd_doc['RXD_CATXRN']))
        references.append(ref)
        
    context['references'] = sorted(references, key=lambda x: x['y'] if x['y'] != 'unk' else -1, reverse=True)
    from makeit.retro.conditions import average_template_list
    context['suggested_conditions'] = average_template_list(INSTANCE_DB, CHEMICAL_DB, reference_ids)

    return render(request, 'template.html', context)

@login_required
def rxid_target(request, rxid):
    '''
    Examines a reaction record from its id in the database
    where rxid == str(rxid)
    '''
    context = {
        'rxid': rxid
    }

    doc = REACTION_DB.find_one({'_id': int(rxid)})
    if not doc:
        context['err'] = 'RXID not found!'
    else:
        context['rxn_smiles'] = doc['RXN_SMILES']
        context['num_instances'] = doc['RX_NVAR']

    reference_ids = ['{}-{}'.format(rxid, i + 1) for i in range(doc['RX_NVAR'])]
    from makeit.retro.conditions import average_template_list
    context['suggested_conditions'] = average_template_list(INSTANCE_DB, CHEMICAL_DB, reference_ids)

    return render(request, 'rxid.html', context)

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
    data = {'err': False}
    smiles = request.GET.get('smiles', None)
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

@login_required
def draw_smiles(request, smiles):
    '''
    Returns a png response for a target smiles
    '''
    from makeit.retro.draw import MolsSmilesToImage
    response = HttpResponse(content_type = 'img/png')
    MolsSmilesToImage(str(smiles)).save(response, 'png')
    return response

@login_required
def draw_smiles_page(request, smiles):
    '''
    Same as draw_smiles but loads as an indepdent page

    OUTDATED: USE PRIMARY /draw/ PAGE (8-31-2016)
    '''
    context = {
        'image_url': reverse('draw_smiles', kwargs={'smiles':smiles}),
        'label_title': 'SMILES',
        'label': smiles,
    }
    if request.method == 'POST':
        context['form'] = DrawingInputForm(request.POST)
    else:
        context['form'] = DrawingInputForm()
    return render(request, 'image.html', context)

@login_required
def draw_template(request, template):
    '''
    Returns a png response for a reaction SMARTS template
    '''
    from makeit.retro.draw import TransformStringToImage
    response = HttpResponse(content_type = 'img/png')
    TransformStringToImage(str(template)).save(response, 'png')
    return response

@login_required
def draw_template_page(request, template):
    '''
    Same as draw_template but loads as an indepdent page

    OUTDATED: USE PRIMARY /draw/ PAGE (8-31-2016)
    '''
    context = {
        'image_url': reverse('draw_template', kwargs={'template':template}),
        'label_title': 'SMARTS template',
        'label': template,
    }
    if request.method == 'POST':
        context['form'] = DrawingInputForm(request.POST)
    else:
        context['form'] = DrawingInputForm()
    return render(request, 'image.html', context)

@login_required
def draw_reaction(request, smiles):
    '''
    Returns a png response for a SMILES reaction string
    '''
    from makeit.retro.draw import ReactionStringToImage
    response = HttpResponse(content_type = 'img/png')
    ReactionStringToImage(str(smiles)).save(response, 'png')
    return response

@login_required
def draw_reaction_page(request, smiles):
    '''
    Same as draw_reaction but loads as an independent page

    OUTDATED: USE PRIMARY /draw/ PAGE (8-31-2016)
    '''
    context = {
        'image_url': reverse('draw_reaction', kwargs={'smiles':smiles}),
        'label_title': 'Reaction SMILES',
        'label': smiles,
    }
    if request.method == 'POST':
        context['form'] = DrawingInputForm(request.POST)
    else:
        context['form'] = DrawingInputForm()
    return render(request, 'image.html', context)

@login_required
def draw(request):
    '''
    Landing page for al draw_*_page functions
    '''
    context = {}

    if request.method == 'POST':
        context['form'] = DrawingInputForm(request.POST)
        if not context['form'].is_valid():
            context['err'] = 'Could not parse!'
        else:
            # Identify target
            text = context['form'].cleaned_data['text']
            try:
                if 'mol' in request.POST:
                    #text = resolve_smiles(text)
                    context['image_url'] = reverse('draw_smiles', kwargs={'smiles':text})
                    context['label_title'] = 'Molecule SMILES'
                    context['label'] = text
                elif 'rxn' in request.POST:
                    #text = '>>'.join([resolve_smiles(frag) for frag in text.split('>>')])
                    context['image_url'] = reverse('draw_reaction', kwargs={'smiles':text})
                    context['label_title'] = 'Reaction SMILES'
                    context['label'] = text
                elif 'tform' in request.POST:
                    context['image_url'] = reverse('draw_template', kwargs={'template':text})
                    context['label_title'] = 'Template SMARTS'
                    context['label'] = text
                else:
                    context['err'] = 'Did not understand request'

            except Exception as e:
                context['err'] = e

    else:
        context['form'] = DrawingInputForm()

    return render(request, 'image.html', context)

@login_required
def draw_synthesis_tree(request, target = None, id = ''):
    '''
    Draw a synthesis tree
    '''

    context = {}

    def assign_ids(target, counter = 1):
        if target['children'] == []:
            target['id'] = counter
            return counter + 1
        for i in range(len(target['children'])):
            counter = assign_ids(target['children'][i], counter)
            counter += 1
        target['id'] = counter
        return counter + 1


    def score_target(target):
        if target['children'] == []:
            if 'is_reaction' in target: raise ValueError('Reaction nodes need children!')
            target['score'] = score_smiles(target['smiles'], ppg = target['ppg'])
        else:
            # Recursively score
            if 'is_reaction' in target: # reactions incur fixed -100 cost
                target['score'] = -100 + sum([score_target(child) for child in target['children']])
            else: # a chemical doesn't add any cost
                target['score'] = sum([score_target(child) for child in target['children']])
        return target['score']

    def score_smiles(smiles, ppg = 0):
        if ppg != 0: return 0 # buyable molecules are free! (ish)
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return 1000 # invalid molecule, infinite score
        total_atoms = mol.GetNumHeavyAtoms()
        ring_bonds = sum([b.IsInRing() - b.GetIsAromatic() for b in mol.GetBonds()])
        chiral_centers = len(Chem.FindMolChiralCenters(mol))
        return  	- 2.00 * float(total_atoms) ** 1.5 \
                    - 1.00 * float(ring_bonds) ** 1.5 \
                    - 0.00 * float(chiral_centers) ** 2.0

    def chem_dict(smiles, children = []):
        return {
            'is_chemical': True,
            'smiles' : smiles,
            'img' : reverse('draw_smiles', kwargs={'smiles':smiles}),
            'ppg' : Pricer.lookup_smiles(smiles),
            'children': children,
        }

    def rxn_dict(info, children = []):
        return {
            'is_reaction': True,
            'info': info,
            'children': children,
        }

    # def id_in_target(target, id):
    # 	if target['id'] == int(id):
    # 		return True
    # 	if target['children'] == []:
    # 		return False
    # 	return any([id_in_target(child, id) for child in target['children']])

    # def expand_by_id(target, id):
    # 	'''Given a target tree and an ID, recursively find it and expand or collapse'''
    # 	if target['id'] == int(id):
    # 		if target['children'] == []:
    # 			precursors = RetroTransformer.perform_retro(smiles).return_top(n=50)
    # 			for smiles in precursor['smiles_split']
    # 			pass
    # 		else:
    # 			#collapse
    # 			pass
    # 	else:
    # 		[collapse_or_expand_by_id(child, id) for child in target['children']]

    # Require a set of trees to use "clicked_on_id"
    if id and 'working_trees' not in request.session:
        context['err'] = 'There is not a set of synthesis trees in memory, so you cannot select an ID!'
        return render(request, 'tree.html', context)

    if 'working_trees' in request.session:
        # Get trees and expand/collapse as necessary
        trees = request.session['working_trees']

        # for i in range(len(trees)):
        # 	if id_in_target(trees[i], id):



    else:
        # Give an example
        target = 'CCCOCCC'
        tree1 = chem_dict(target, children = [
            rxn_dict('rxn1', children = [
                chem_dict('CCCO'),
                chem_dict('CCC[Br]')
            ]),
        ])
        tree2 = chem_dict(target, children = [
            rxn_dict('rxn2', children = [
                chem_dict('CCCO'),
                chem_dict('CCC[Cl]', children = [
                    rxn_dict('Needs source of [Cl]', children = [
                        chem_dict('CCCO'),
                    ])
                ])
            ]),
        ])
        target += ' (example)'
        trees = [tree2, tree1]

    counter = 1
    for tree in trees:
        counter = assign_ids(tree, counter)
        score_target(tree)
    context['target'] = trees[0]['smiles']
    context['trees'] = sorted(trees, key = lambda x: x['score'], reverse = True)[:50]

    request.session['working_trees'] = context['trees']

    return render(request, 'tree.html', context)

def nn_predictor_setup(request,):
    context = {}
    if request.method == 'POST':
        context['form'] = nnSetup(request.POST)

        if not context['form'].is_valid():
            context['err'] = 'Invalid input'
        else:
            nnPred = NNConditionPredictor()
            nnPred.load_predictor(context['form'].cleaned_data)
            ### **** reactions, paths ****

    else:
        context['form'] = nnSetup()
    return render(request, 'nnSetupInput.html', context)

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
def draw_fig(request, fig):
    '''
    Returns a png response for figure object
    '''
    response = HttpResponse(content_type = 'img/png')
    canvas = FigureCanvas(fig)
    canvas.print_png(response)
    return response


def sep_input(request):
    context = {}
    if request.method == 'POST':
        context['form'] = sepInput(request.POST)

        if not context['form'].is_valid():
            context['err'] = 'Invalid input'
        else:
            SD = SeparationDesigner()
            SD.load_designer(context['form'].cleaned_data)
            output = SD.opt_sep_pH_flow()

            # Convert Figure objects into URLs
            imgs = []
            for f in output[1:]:
                imgs.append(reverse('draw_fig', kwargs={'fig': f}))

            # Pass output into context
            context['output'] = {
                'result': output[0],
                'imgs': imgs
            }
    else:
        context['form'] = sepInput()
    return render(request, 'sepInput.html', context)
