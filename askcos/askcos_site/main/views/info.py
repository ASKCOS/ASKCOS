from django.shortcuts import render, HttpResponse, redirect
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.conf import settings
import django.contrib.auth.views
from pymongo.message import bson
from bson.objectid import ObjectId
from pymongo import MongoClient
from pymongo.errors import ServerSelectionTimeoutError
import makeit.global_config as gc

from .users import can_view_reaxys

client = MongoClient(
    gc.MONGO['path'],
    gc.MONGO['id'],
    connect=gc.MONGO['connect']
)
db_name = gc.RETRO_TEMPLATES['database']
collection = gc.RETRO_TEMPLATES['collection']
retro_templates = client[db_name][collection]

def template_target_export(request, id):
    rx_ids = template_target(request, id, return_refs_only=True)

    txt = '{"fileName":"reaxys_query.json", "version":"1.0","content": {"id":"root","facts": ['
    for i, rx_id in enumerate(rx_ids):
        if i != 0:
            txt += ','
        txt += '{"id":"Reaxys487",'
        if i != 0:
            txt += '"logicOperator":"OR",'
        txt += '"fields": [{"value":"%i","boundOperator":"op_num_equal","id":"RX.ID","displayName":"Reaction ID"}],' % (rx_id)
        txt += '"fieldsLogicOperator":"AND","exist":false,"bio":false}'

    txt += '], "exist":false,"bio":false } }'

    response = HttpResponse(txt, content_type='text/csv')
    response['Content-Disposition'] = 'attachment;filename=reaxys_query.json'
    return response

#@login_required
def template_target(request, id, return_refs_only=False):
    '''
    Examines a template from its id in the database
    where id == str(_id)

    Reaxys templates refer to instances, but we need the
    reactions for visualization

    If return_refs_only, builds .json query for Reaxys
    '''
    context = {}

    transform = retro_templates.find_one({'_id': ObjectId(id)})

    if not transform:
        context['err'] = 'Transform not found'
        return render(request, 'template.html', context)

    context['template'] = transform
    context['template']['id'] = id
    reference_ids = transform['references']

    if return_refs_only:
        return list(sorted(set(int(_id.split('-')[0]) for _id in reference_ids)))[:500]

    references = []
    rx_docs = {}; xrn_to_smiles = {}
    context['total_references'] = len(reference_ids)

    context['cannot_view_any'] = True
    context['ref_ids'] = '; '.join([str(y) for y in sorted(set(int(_id.split('-')[0]) for _id in reference_ids))])

    return render(request, 'template.html', context)
