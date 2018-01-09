from django.shortcuts import render, HttpResponse, redirect
from django.template.loader import render_to_string
from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.conf import settings
import django.contrib.auth.views
import json
import os

import rdkit.Chem as Chem 
import urllib2

from ...askcos_celery.contextrecommender.worker import get_context_recommendation

from ..utils import get_name_from_smiles
from .users import can_control_robot

@login_required
def export_synth_results(request):
    if 'last_synth_interactive' not in request.session:
        return index(request, err='Could not find synth results to save?')
    
    synth_data = request.session['last_synth_interactive']
    txt = '%s\r\n' % SYNTH_FOOTNOTE
    txt += 'Reactants: %s\r\n' % synth_data['reactants']
    txt += 'Reagents: %s\r\n' % synth_data['reagents']
    txt += 'Temperature: %s\r\n' % synth_data['temperature']
    txt += 'Solvent: %s\r\n' % synth_data['solvent']
    txt += '\r\n'
    txt += '%s\t%s\t%s\t%s\r\n' % ('Rank', 'SMILES', 'Probability', 'Score')
    for outcome in synth_data['outcomes']:
        txt += '%s\t%s\t%s\t%s\r\n' % (outcome['rank'], outcome['smiles'], outcome['prob'], outcome['score'])
    response = HttpResponse(txt, content_type='text/csv')
    response['Content-Disposition'] = 'attachment;filename=export.csv'
    return response
    
@login_required
def export_retro_results(request, _id=1):
    if 'last_retro_interactive' not in request.session:
        return index(request, err='Could not find retro results to save?')
    if not can_control_robot(request):
        return index(request, err='You do not have permission to be here!')
    
    try:
        _id = int(_id)
    except Exception:
        return index(request, err='Passed _id was non-integer...?')

    if _id > len(request.session['last_retro_interactive']):
        return index(request, err='Requested index out of range?')

    tree = request.session['last_retro_interactive'][_id - 1]

    class mutable:
        bays = []; 
        reagents = []
        bay_id = 1; 
        reagent_id = 1 # start at 1 and work up, then change at the end

    def recursive_get_bays(chem):
        # Is this the product of a reaction? Then start defining the reactor
        if chem['children']:
            this_bay = {'type': 'reactor', 'id': mutable.bay_id}
            mutable.bay_id += 1
            
            # Outlet is this chemical (smiles)
            this_bay['outlet_smiles'] = chem['smiles']
            this_bay['outlet'] = get_name_from_smiles(chem['smiles'])
            # Only child is a reaction
            rxn = chem['children'][0]
            # Make sure this is not a convergent synthesis
            if sum([len(child['children']) > 0 for child in rxn['children']]) > 1:
                return index(request, err='Cannot perform convergent syntheses yet...')

            # Save reaction SMILES for the reactor
            this_bay['reaction_smiles'] = rxn['smiles']
            this_bay['necessary_reagent'] = rxn['necessary_reagent']
            
            # Get context recommendations (just one for now, no forward evaluator)
            res = get_context_recommendation.delay(rxn['smiles'], n=10)
            contexts = res.get(60)
            if not contexts:
                # use default conditions
                (T1, slvt1, rgt1, cat1) = 20., 'C1COCC1', '', ''
                slvt1_name = 'THF (unk. solvent)'
                this_bay['context_failed'] = True
            else:
                (T1, slvt1, rgt1, cat1, t1, y1) = contexts[0]
                slvt1_name = get_name_from_smiles(slvt1)
                this_bay['context_failed'] = False
            this_bay['temperature'] = T1
            this_bay['reaction_solvent_smiles'] = slvt1
            this_bay['reaction_solvent_name'] = slvt1_name

            # Define inlets
            this_bay['inlets'] = []

            # Define inlets - reagents/catalysts
            for rgt in rgt1.split('.') + cat1.split('.'):
                if not rgt:
                    continue
                reagent = {
                    'id': mutable.reagent_id,
                    'smiles': rgt,
                    'name': '%s in %s' % (get_name_from_smiles(rgt), slvt1_name),
                    'type': 'reagent',
                }
                mutable.reagent_id += 1
                mutable.reagents.append(reagent)
                this_bay['inlets'].append(reagent)

            # Define inlets - new reactants/feeds
            for child in rxn['children']:
                if not child['children']:
                    reagent = {
                        'id': mutable.reagent_id,
                        'smiles': child['smiles'],
                        'name': '%s in %s' % (get_name_from_smiles(child['smiles']), slvt1_name),
                        'type': 'reactant',
                        'ppg': child['ppg'],
                    }
                    mutable.reagent_id += 1
                    mutable.reagents.append(reagent)
                    this_bay['inlets'].append(reagent)

            # Save this bay
            mutable.bays.append(this_bay)

            # Now add the next bays as needed
            for child in rxn['children']:
                if child['children']:
                    recursive_get_bays(child)

    # Actually perform this on the current head node
    recursive_get_bays(tree)
    bays = mutable.bays 
    reagents = mutable.reagents

    # Renumber bays
    num_bays = len(bays)
    for bay in bays:
        bay['id'] = num_bays - bay['id'] + 1
    # Renumber reagents
    num_reagents = len(reagents)
    for reagent in reagents:
        reagent['id'] = num_reagents - reagent['id'] + 1

    # Report (for debugging mostly)
    print('{} bays'.format(len(bays)))
    for bay in sorted(bays, key=lambda x:int(x['id'])):
        print('\nBay ID {:2d}'.format(bay['id']))
        print('    {:20s}{:20s}'.format('type', bay['type']))
        print('    {:20s}{:20s}'.format('solvent', bay['reaction_solvent_name']))
        print('    {:20s}{:20s}'.format('temperature', str(bay['temperature'])))
        inlets = ', '.join([str(rgt['id']) for rgt in bay['inlets']])
        if bay['id'] != 1:
            inlets += ', and bay {} outlet'.format(bay['id'] - 1)
        print('    {:20s}{:20s}'.format('inlets', inlets))
        print('    {:20s}{:20s}'.format('outlet', bay['outlet']))
    print('\n{} reagents'.format(len(reagents)))
    for reagent in sorted(reagents, key=lambda x:int(x['id'])):
        print('ID {:2d}: {}'.format(reagent['id'], reagent['name']))

    

    # Reformat for Suzanne
    bay_dict = {}
    for bay in bays:
        key = 'bay{}'.format(bay['id'])
        bay_dict[key] = {
            'transformation': bay['reaction_smiles'],
            'temp': bay['temperature'],
        }
        for i in range(len(bay['inlets'])):
            bay_dict[key]['reagent{}'.format(i+1)] = bay['inlets'][i]['name']
    
    bays = json.dumps(bays)
    print('Post-json dump: {}'.format(bays))

    bay_dict = json.dumps(bay_dict)
    print('Bay-dict for Suzanne: {}'.format(bay_dict))


    return post_data_to_external_page(request, 'http://Apex-env.qr9fgjkxpf.us-east-2.elasticbeanstalk.com/move_robot/', 
        {'bay_dict': bay_dict, 'bays': bays})

@login_required
def post_data_to_external_page(request, url, data={}):
    '''
    Returns a dummy page that immediately posts to an external site
    '''
    return render(request, 'redirect.html', {'data': data, 'url': url})
