from bson.objectid import ObjectId
from django.http import JsonResponse
import makeit.global_config as gc
from pymongo import MongoClient

client = MongoClient(
    gc.MONGO['path'],
    gc.MONGO['id'],
    connect=gc.MONGO['connect']
)
db_name = gc.RETRO_TRANSFORMS_CHIRAL['database']
collection = gc.RETRO_TRANSFORMS_CHIRAL['collection']
retro_chiral = client[db_name][collection]

def template(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    _id = request.GET.get('id')
    transform = retro_chiral.find_one({'_id': ObjectId(_id)})
    if not transform:
        transform = retro_achiral.find_one({'_id': ObjectId(_id)})
    if not transform:
        transform = synth.find_one({'_id': ObjectId(_id)})
    if not transform:
        resp['error'] = 'Cannot find template with id '+_id
        return JsonResponse(resp)
    transform['_id'] = _id
    transform.pop('product_smiles', None)
    transform.pop('name', None)
    refs = transform.pop('references', [''])
    transform['references'] = list(map(lambda x: x.split('-')[0], refs))
    resp['template'] = transform
    return JsonResponse(resp)
