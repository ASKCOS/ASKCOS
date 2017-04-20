from django.shortcuts import render, HttpResponse, redirect
from django.template.loader import render_to_string
from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.conf import settings
import django.contrib.auth.views
from forms import SmilesInputForm, DrawingInputForm, is_valid_smiles
from bson.objectid import ObjectId
import time

import rdkit.Chem as Chem 
import urllib2

from askcos_site.main.globals import RetroTransformer, RETRO_FOOTNOTE, SynthTransformer, SYNTH_FOOTNOTE, REACTION_DB, INSTANCE_DB, CHEMICAL_DB, BUYABLE_DB, SOLVENT_DB, Pricer, TransformerOnlyKnown, builder, predictor, PREDICTOR_FOOTNOTE
from forms import nnSetup, sepInput
from askcos_site.functions.SPARC import SeparationDesigner
from askcos_site.functions.nnPredictor import NNConditionPredictor

def the_time(request):
    context = {
        'the_time': time.time(),
    }

    return render(request, 'get_time.html', context)

def get_the_time(request):
    prevTime = request.GET.get('prevTime', None)
    print('Prev time: {}'.format(prevTime))
    data = {
        'newTime': time.time()
    }
    return JsonResponse(data)

def index(request):
    '''
    Homepage
    '''
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
            if 'retro_lit' in request.POST: return redirect('retro_lit_target', smiles = smiles)
            if 'retro' in request.POST:
                if is_valid_smiles(smiles):
                    return redirect('retro_target', smiles = smiles)
                else:
                    context['err'] = 'Invalid SMILES string: {}'.format(smiles)
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
        {'name': 'Itraconazole', 'smiles': 'CCC(C)N1N=CN(C1=O)c2ccc(cc2)N3CCN(CC3)c4ccc(OC[C@H]5CO[C@@](Cn6cncn6)(O5)c7ccc(Cl)cc7Cl)cc4'}
    ]

    context['footnote'] = RETRO_FOOTNOTE
    return render(request, 'retro.html', context)

@login_required
def retro_target(request, smiles, max_n = 50):
    '''
    Given a target molecule, render page
    '''

    # Render form with target
    context = {}
    context['form'] = SmilesInputForm({'smiles': smiles})

    # Look up target
    smiles_img = reverse('draw_smiles', kwargs={'smiles':smiles})
    context['target'] = {
        'smiles': smiles,
        'img': smiles_img
    }

    # Perform retrosynthesis
    startTime = time.time()
    result = RetroTransformer.perform_retro(smiles)
    context['precursors'] = result.return_top(n = 50)
    elapsedTime = time.time() - startTime

    # Change 'tform' field to be reaction SMARTS, not ObjectID from Mongo
    # Also add up total number of examples
    for (i, precursor) in enumerate(context['precursors']):
        context['precursors'][i]['tforms'] = \
            [dict(RetroTransformer.lookup_id(_id), **{'id':str(_id)}) for _id in precursor['tforms']]
        context['precursors'][i]['mols'] = []
        for smiles in precursor['smiles_split']:
            ppg = Pricer.lookup_smiles(smiles, alreadyCanonical = True)
            context['precursors'][i]['mols'].append({
                'smiles': smiles,
                'ppg': '${}/g'.format(ppg) if ppg else 'cannot buy'
            })

    context['footnote'] = RETRO_FOOTNOTE + ', generated results in %0.3f seconds' % elapsedTime
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
    context['warn'] = 'This module currently uses a shared backend object in Django, so there will be conflicts if more than one person tries to use this page at the same time.'

    return render(request, 'retro_interactive.html', context)

@login_required
def synth_interactive(request):
    '''Builds an interactive forward synthesis page'''

    context = {}
    context['warn'] = 'This module currently uses a shared backend object in Django, so there will be conflicts if more than one person tries to use this page at the same time.'
    context['footnote'] = PREDICTOR_FOOTNOTE

    solvent_choices = []
    for doc in SOLVENT_DB.find({'_id': {'$ne': 'default'}}):
        solvent_choices.append({
            'smiles': doc['smiles'],
            'name': doc['name'],
        })
    context['solvent_choices'] = sorted(solvent_choices, key = lambda x: x['name'])
    context['mincount_default'] = settings.SYNTH_TRANSFORMS['mincount']

    return render(request, 'synth_interactive.html', context)

def ajax_smiles_to_image_retro(request):
    '''Takes an Ajax call with a smiles string
    and returns the HTML for embedding an image'''

    smiles = request.GET.get('smiles', None)
    print('SMILES from Ajax: {}'.format(smiles))

    # Make sure we are not still building a previous tree...
    builder.stop_building(timeout = 3) # just to be sure

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        # Try to resolve using NIH
        try:
            smiles = urllib2.urlopen('https://cactus.nci.nih.gov/chemical/structure/{}/smiles'.format(smiles)).read()
            print('Resolved smiles -> {}'.format(smiles))
            mol = Chem.MolFromSmiles(smiles)
            if not mol: return JsonResponse({'err': True})
        except:
            return JsonResponse({'err': True})

    url = reverse('draw_smiles', kwargs={'smiles':smiles})
    data = {
        'html': '<img src="' + url + '">',
        'smiles': smiles,
    }

    return JsonResponse(data)

def ajax_smiles_to_image_synth(request):
    '''Takes an Ajax call with a smiles string
    and returns the HTML for embedding an image'''
    data = {'err': False}

    smiles = request.GET.get('smiles', None)
    print('SMILES from Ajax: {}'.format(smiles))

    if not smiles:
        data['html'] = 'no reagents'
        data['smiles'] = smiles
        return JsonResponse(data)

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        # Try to resolve using NIH
        try:
            smiles = smiles.replace('.', ' and ')
            smiles_list = [x.strip() for x in smiles.split(' and ')]
            for i, smi in enumerate(smiles_list):
                print(smi)
                smi = urllib2.urlopen('https://cactus.nci.nih.gov/chemical/structure/{}/smiles'.format(urllib2.quote(smi))).read()
                print('Resolved smiles -> {}'.format(smi))
                mol = Chem.MolFromSmiles(smi)
                if not mol: return JsonResponse({'err': True})
                smiles_list[i] = smi
            smiles = '.'.join(smiles_list)
        except:
            return JsonResponse({'err': True})

    print('Resolved smiles: {}'.format(smiles))

    url = reverse('draw_smiles', kwargs={'smiles':smiles})
    data = {
        'html': '<img src="' + url + '">',
        'smiles': smiles,
    }
    return JsonResponse(data)


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
    predictor.set_context(T = temperature, reagents = reagents, solvent = solvent)
    print('Set context')
    predictor.run_in_foreground(reactants = reactants, intended_product = '', quit_if_unplausible = False, mincount = mincount)
    print('Ran in foreground')
    outcomes = predictor.return_top(n = 25)
    print('Got top outcomes, length {}'.format(len(outcomes)))
    data['html_time'] = '{:.3f} seconds elapsed'.format(time.time() - startTime)

    if outcomes:
        data['html'] = render_to_string('synth_outcomes_only.html', {'outcomes': outcomes})
    else:
        data['html'] = 'No outcomes found? That is weird...'
    return JsonResponse(data)


def ajax_start_retro(request):
    '''Start builder'''
    smiles = request.GET.get('smiles', None)
    max_depth = int(request.GET.get('max_depth', 4))
    max_branching = int(request.GET.get('max_branching', 25))
    data = {'err': False}
    if builder.is_running() and builder.is_target(smiles):
        builder.unpause()
    else:
        builder.stop_building(timeout = 3) # just to be sure
        builder.start_building(smiles, max_depth = max_depth, max_branching = max_branching)
    return JsonResponse(data)

def ajax_pause_retro(request):
    '''Pause builder'''
    smiles = request.GET.get('smiles', None)
    data = {'err': False}
    if builder.is_target(smiles):
        builder.pause()
    else:
        data['err'] = True
        data['message'] = 'Cannot pause if we have not started running'
    return JsonResponse(data)

def ajax_stop_retro(request):
    '''Stop builder'''
    data = {'err': False}
    if builder.is_running():
        builder.stop_building(timeout = 3)
    else:
        data['err'] = True
        data['message'] = 'Cannot stop if we arent running!'
    return JsonResponse(data)

def ajax_update_retro_stats(request):
    '''Update the statistics only'''
    smiles = request.GET.get('smiles', None)
    data = {'err': False}
    if builder.is_target(smiles):
        data['html_stats'] = builder.info_string()
    else:
        data['err'] = True
        data['message'] = 'Cannot update stats if we have not started running!'
    return JsonResponse(data)

def ajax_update_retro(request):
    '''Update displayed results'''
    data = {'err': False}
    smiles = request.GET.get('smiles', None)
    if builder.is_target(smiles):
        data['html_stats'] = builder.info_string()
        print(builder.info_string())
        trees = builder.get_trees_iddfs()
        print('Got trees')
        if trees:
            data['html_trees'] = render_to_string('trees_only.html', {'trees': trees})
        else:
            data['html_trees'] = 'No trees resulting in buyable chemicals found!'
    else:
        data['err'] = True
        data['message'] = 'Cannot show results if we have not started running'
    return JsonResponse(data)

@login_required
def synth(request):
    '''
    Forward synthesis homepage
    '''
    context = {}

    if request.method == 'POST':
        context['form'] = SmilesInputForm(request.POST)
        if not context['form'].is_valid():
            context['err'] = 'Could not parse!'
        else:
            # Identify target
            smiles = context['form'].cleaned_data['smiles']
            if is_valid_smiles(smiles):
                return redirect('synth_target', smiles = smiles)
            else:
                context['err'] = 'Invalid SMILES string: {}'.format(smiles)
    else:
        context['form'] = SmilesInputForm()

    context['footnote'] = SYNTH_FOOTNOTE
    return render(request, 'synth.html', context)

@login_required
def synth_target(request, smiles, max_n = 50):
    '''
    Given a set of reactants as a single SMILES string, find products
    '''

    # Render form with target
    context = {}
    context['form'] = SmilesInputForm({'smiles': smiles})

    # Look up target
    smiles_img = reverse('draw_smiles', kwargs={'smiles':smiles})
    context['target'] = {
        'smiles': smiles,
        'img': smiles_img,
    }

    # Perform forward synthesis
    result = SynthTransformer.perform_forward(smiles)
    context['products'] = result.return_top(n = 50)
    #print(context['products'])
    # Change 'tform' field to be reaction SMARTS, not ObjectID from Mongo
    # (add "id" field to be string version of ObjectId "_id")
    for (i, product) in enumerate(context['products']):
        context['products'][i]['tforms'] = \
            [dict(SynthTransformer.lookup_id(_id), **{'id':str(_id)}) for _id in product['tforms']]

    context['footnote'] = SYNTH_FOOTNOTE
    return render(request, 'synth.html', context)

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

    references = []; docs = []
    context['total_references'] = len(reference_ids)
    for i, ref_id in enumerate(reference_ids):
        doc = None
        try:
            doc = REACTION_DB.find_one({'_id': int(ref_id.split('-')[0])})
        except Exception as e:
            print(e)
        if doc:
            docs.append(doc)
            references.append({
                'label': ref_id,
                'reaction_smiles': doc['RXN_SMILES'],
            })
        else:
            print('Could not find doc with example')

    from makeit.retro.conditions import average_template_list
    context['suggested_conditions'] = average_template_list(INSTANCE_DB, CHEMICAL_DB, reference_ids)

    context['references'] = references[:15]
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
    ppg = Pricer.lookup_smiles(smiles, alreadyCanonical = True)
    if ppg:
        response.write('${}/g'.format(ppg))
    else:
        response.write('cannot buy')
    return response

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
                    context['image_url'] = reverse('draw_smiles', kwargs={'smiles':text})
                    context['label_title'] = 'Molecule SMILES'
                    context['label'] = text
                elif 'rxn' in request.POST:
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