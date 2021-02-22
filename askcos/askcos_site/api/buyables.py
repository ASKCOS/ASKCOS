import json
from django.http import JsonResponse
from django.contrib.auth.decorators import login_required
from makeit import global_config as gc
from askcos_site.main.views.users import can_modify_buyables
from pymongo import MongoClient
from bson import ObjectId
import pandas as pd
import io
from rdkit import Chem

mongo = MongoClient(gc.MONGO['path'])
buyables_db = mongo[gc.BUYABLES['database']][gc.BUYABLES['collection']]


def add_buyable_list_to_db(buyable_list, allow_overwrite=True):
    result = {
        'error': None,
        'count': 0,
        'added': [],
        'updated': [],
        'added_count': 0,
        'updated_count': 0,
        'duplicate_count': 0,
        'total': len(buyable_list)
    }
    for buyable in buyable_list:
        res = add_buyable_to_db(buyable, allow_overwrite=allow_overwrite)
        if not res['error']:
            if res.get('duplicate'):
                result['duplicate_count'] += 1
            if res.get('inserted'):
                result['added'].append(res['inserted'])
                result['added_count'] += 1
            if res.get('updated'):
                result['updated'].append(res['updated'])
                result['updated_count'] += 1
            result['count'] += 1
    return result

def add_buyable_to_db(buyable, allow_overwrite=True):
    result = {'error': None}
    smiles = buyable.get('smiles')
    ppg = float(buyable.get('ppg', 0.0))
    source = buyable.get('source')
    if not smiles:
        result['error'] = 'Smiles not provided for buyables entry'
        return result
    if not ppg:
        result['error'] = 'Price not provided. Price must be a nonzero float value.'
        return result
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        result['error'] = 'Cannot parse smiles with rdkit'
        return result
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    new_doc = {
        'smiles': smiles,
        'ppg': ppg,
        'source': source
    }
    existing_doc = buyables_db.find_one({'smiles': smiles})
    if existing_doc and allow_overwrite:
        buyables_db.update_one(
            {'smiles': smiles},
            {'$set': {'ppg': ppg, 'source': source}}
        )
        new_doc['_id'] = str(existing_doc['_id'])
        result['updated'] = new_doc
    elif existing_doc and not allow_overwrite:
        result['duplicate'] = True
    else:
        insert_result = buyables_db.insert_one(new_doc)
        if not insert_result.inserted_id:
            result['error'] = 'Addition of buyable failed'
            return result
        new_doc['_id'] = str(new_doc['_id'])
        result['inserted'] = new_doc
    return result


def buyables(request):
    search = request.GET.get('q', None)
    source = request.GET.get('source', None)
    regex = request.GET.get('regex', None) in ['true', 'True']
    limit = int(request.GET.get('limit', 100))
    canon = request.GET.get('canonicalize', 'True') in ['True', 'true']
    query = {}
    if search:
        if regex:
            query['smiles'] = {'$regex': '.*{}.*'.format(search)}
        else:
            if canon:
                mol = Chem.MolFromSmiles(search)
                if mol:
                    search = Chem.MolToSmiles(mol, isomericSmiles=True)
            query['smiles'] = search
    if source:
        query['source'] = source
    resp = {}
    search_result = list(buyables_db.find(
        query,
        {'smiles': 1, 'ppg': 1, 'source': 1}
    ).limit(limit))
    for doc in search_result:
        doc['_id'] = str(doc['_id'])
    resp['buyables'] = search_result
    resp['search'] = search
    return JsonResponse(resp)

def delete_buyable(request):
    resp = {}
    if not can_modify_buyables(request):
        resp['error'] = 'You are not authorized to modify the buyables database. Please ask your site administrator for permission.'
        return JsonResponse(resp)
    buyable_id = request.GET.get('_id', None)
    delete_result = buyables_db.delete_one({'_id': ObjectId(buyable_id)})
    if delete_result.deleted_count != 1:
        resp['error'] = 'Could not find buyable'
        return JsonResponse(resp)
    resp['success'] = True
    return JsonResponse(resp)

def upload_buyable(request):
    resp = {}
    if not can_modify_buyables(request):
        resp['error'] = 'You are not authorized to modify the buyables database. Please ask your site administrator for permission.'
        return JsonResponse(resp)
    print(request.FILES)
    file_format = request.POST.get('format')
    upload_file = request.FILES.get('file')
    return_limit = int(request.POST.get('returnLimit', 1000))
    allow_overwrite = request.POST.get('allowOverwrite', 'True') in ['True', 'true']
    try:
        content = upload_file.read()
    except:
        resp['error'] = 'Cannot read file!'
        return JsonResponse(resp)
    if file_format == 'json':
        try:
            upload_json = json.loads(content.decode('utf-8'))
        except json.JSONDecodeError:
            resp['error'] = 'Cannot parse json!'
            return JsonResponse(resp)
    elif file_format == 'csv':
        try:
            df = pd.read_csv(io.StringIO(content.decode('utf-8')))
            upload_json = df.to_dict(orient='records')
        except Exception as e:
            print(e)
            resp['error'] = 'Cannot parse csv!'
            return JsonResponse(resp)
    else:
        resp['error'] = 'File format not supported!'
        return JsonResponse(resp)
    if not isinstance(upload_json, list) or len(upload_json) == 0 or not isinstance(upload_json[0], dict):
        resp['error'] = 'Improperly formatted json!'
        return JsonResponse(resp)
    result = add_buyable_list_to_db(upload_json, allow_overwrite=allow_overwrite)
    if result.get('error'):
        resp['error'] = result['error']
    resp['added'] = result['added'][:return_limit]
    resp['updated'] = result['updated'][:return_limit]
    resp['added_count'] = result['added_count']
    resp['updated_count'] = result['updated_count']
    resp['duplicate_count'] = result['duplicate_count']
    resp['count'] = result['count']
    resp['total'] = result['total']
    return JsonResponse(resp)

def add_buyable(request):
    resp = {}
    if not can_modify_buyables(request):
        resp['error'] = 'You are not authorized to modify the buyables database. Please ask your site administrator for permission.'
        return JsonResponse(resp)
    smiles = request.GET.get('smiles', None)
    ppg = float(request.GET.get('ppg', 0.0))
    source = request.GET.get('source', '')
    allow_overwrite = request.GET.get('allowOverwrite', 'True') in ['True', 'true']
    result = add_buyable_to_db({
        'smiles': smiles,
        'ppg': ppg,
        'source': source
    }, allow_overwrite=allow_overwrite)
    if result['error']:
        resp['error'] = result['error']
        return JsonResponse(resp)
    resp['success'] = True
    if result.get('duplicate') == True:
        resp['duplicate'] = True
    if result.get('inserted'):
        resp['inserted'] = result['inserted']
    if result.get('updated'):
        resp['updated'] = result['updated']
    return JsonResponse(resp)