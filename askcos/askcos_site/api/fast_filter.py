from django.http import JsonResponse
from makeit import global_config as gc
from askcos_site.askcos_celery.treebuilder.tb_c_worker import fast_filter_check

def fast_filter(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    reactants = request.GET.get('reactants')
    products = request.GET.get('products')
    res = fast_filter_check.delay(reactants, products)
    outcome = res.get(60)
    score = outcome[0][0]['score']
    resp['score'] = score
    return JsonResponse(resp)