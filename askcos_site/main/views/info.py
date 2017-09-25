from django.shortcuts import render, HttpResponse, redirect
from django.contrib.auth.decorators import login_required
from django.conf import settings
import django.contrib.auth.views
from pymongo.message import bson
from bson.objectid import ObjectId

from ..globals import RetroTransformer, SynthTransformer, \
    REACTION_DB, INSTANCE_DB, CHEMICAL_DB, TEMPLATE_BACKUPS, REACTION_DB_OLD, \
    INSTANCE_DB_OLD, CHEMICAL_DB_OLD

from ..utils import fancyjoin, xrn_lst_to_name_lst

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
    if not transform: 
        transform = SynthTransformer.lookup_id(ObjectId(id))

    # Look through backups...
    if not transform:
        for coll in TEMPLATE_BACKUPS:
            transform = coll.find_one({'_id': ObjectId(id)})
            if transform: 
                context['warn'] = 'This template is out of date (from %s)' % coll._Collection__full_name
                break
        if not transform:
            context['err'] = 'Transform not found, even after looking in backup DBs %s' % ', '.join([coll._Collection__full_name for coll in TEMPLATE_BACKUPS])
            return render(request, 'template.html', context)
    
    context['template'] = transform
    context['template']['id'] = id
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
            'ded': rxd_doc.get('RXD_DED', [''])[0],
        }

        def rxn_lst_to_name_lst(xrn_lst):
            lst = []
            for xrn in xrn_lst:
                if xrn not in xrn_to_smiles: 
                    chem_doc = CHEMICAL_DB.find_one({'_id': xrn})
                    if chem_doc is None:
                        xrn_to_smiles[xrn] = 'Chem-%i' % xrn
                    elif 'IDE_CN' not in chem_doc:
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
        
    context['references'] = sorted(references, key=lambda x: x['y'] if x['y'] != 'unk' else -1, reverse=True)[:500]
    #from makeit.retro.conditions import average_template_list
    #context['suggested_conditions'] = average_template_list(INSTANCE_DB, CHEMICAL_DB, reference_ids)

    return render(request, 'template.html', context)

    
@login_required
def rxid_target(request, rxid):
    '''
    Examines a reaction record from its id in the database
    where rxid == str(rxid)

    UNUSED?
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
