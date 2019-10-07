from django.http import JsonResponse
from django.contrib.auth.decorators import login_required
from celery.result import AsyncResult
from askcos_site.celery import app
from askcos_site.main.models import SavedResults
from pymongo import MongoClient
from makeit import global_config as gc
from celery.result import AsyncResult

client = MongoClient(
    gc.MONGO['path'],
    gc.MONGO['id'],
    connect=gc.MONGO['connect']
)
results_db = client['results']
results_collection = results_db['results']

@login_required
def poll_result(request):
    resp = {}
    id_ = request.GET.get('id')
    resp['id'] = id_
    try:
        result = SavedResults.objects.get(user=request.user, result_id=id_)
    except:
        resp['state'] = -1
    if result:
        resp['state'] = result.result_state
    else:
        resp['state'] = -1
    return JsonResponse(resp)

@login_required
def get_result(request):
    resp = {}
    id_ = request.GET.get('id')
    try:
        result = SavedResults.objects.get(user=request.user, result_id=id_)
    except:
        resp['success'] = False
        resp['error'] = 'No result found!'
        result = None
    if result and result.result_state == 'completed':
        result_doc = results_collection.find_one({'_id': id_})
        resp['result'] = result_doc
    return JsonResponse(resp)

@login_required
def my_results(request):
    resp = {}
    resp['results'] = []
    results = SavedResults.objects.filter(user=request.user)
    for result in results:
        resp['results'].append({
            'id': result.result_id,
            'state': result.result_state,
            'description': result.description,
            'created': result.created.strftime('%B %d, %Y %H:%M:%S %p'),
            'type': result.result_type
        })
    return JsonResponse(resp)

@login_required
def remove_result(request):
    id_ = request.GET.get('id')
    resp = {}
    try:
        result = SavedResults.objects.get(user=request.user, result_id=id_)
    except:
        resp['status'] = -1
    try:
        if result:
            result.delete()
            results_collection.delete_one({'_id': id_})
    except:
        resp['status'] = 0
    resp['status'] = 1
    return JsonResponse(resp)
