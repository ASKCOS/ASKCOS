from rdkit import Chem
from django.http import JsonResponse
from celery.exceptions import TimeoutError
from makeit import global_config as gc
from askcos_site.askcos_celery.treebuilder.tb_c_worker import get_top_precursors as get_top_precursors_c

TIMEOUT = 60
TIMEOUT_ALL = 120

def singlestep(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    run_async = request.GET.get('async', False)
    target = request.GET.get('target')
    template_prioritization = request.GET.get('template_prioritization', gc.relevance)
    precursor_prioritization = request.GET.get('precursor_prioritization', gc.relevanceheuristic)
    mincount = int(request.GET.get('mincount', 0))
    num_templates = int(request.GET.get('num_templates', 100))
    max_cum_prob = float(request.GET.get('max_cum_prob', 0.995))
    apply_fast_filter = request.GET.get('apply_fast_filter', 'True') in ['True', 'true']
    filter_threshold = float(request.GET.get('filter_threshold', 0.75))
    max_branching = int(request.GET.get('num_results', 100))

    if not target:
        resp['error'] = 'Required parameter "target" missing'
        return JsonResponse(resp, status=400)

    mol = Chem.MolFromSmiles(target)
    if not mol:
        resp['error'] = 'Cannot parse target smiles with rdkit'
        return JsonResponse(resp, status=400)
    
    cluster = request.GET.get('cluster', 'True') in ['True', 'true']
    cluster_method = request.GET.get('cluster_method', 'kmeans')
    cluster_feature = request.GET.get('cluster_feature', 'original')
    cluster_fp_type = request.GET.get('cluster_fp_type', 'morgan')
    cluster_fp_length = int(request.GET.get('cluster_fp_length', 512))
    cluster_fp_radius = int(request.GET.get('cluster_fp_radius', 1))

    res = get_top_precursors_c.delay(
        target,
        template_prioritization,
        precursor_prioritization,
        mincount=mincount,
        apply_fast_filter=apply_fast_filter,
        filter_threshold=filter_threshold,
        max_branching=max_branching,
        max_cum_prob=max_cum_prob,
        cluster=cluster,
        cluster_method=cluster_method,
        cluster_feature=cluster_feature,
        cluster_fp_type=cluster_fp_type,
        cluster_fp_length=cluster_fp_length,
        cluster_fp_radius=cluster_fp_radius
    )

    if run_async:
        resp['id'] = res.id
        resp['state'] = res.state
        return JsonResponse(resp)

    if template_prioritization == gc.relevance:
        timeout = TIMEOUT
    else:
        timeout = TIMEOUT_ALL

    try:
        (smiles, precursors) = res.get(timeout)
    except TimeoutError:
        resp['error'] = 'API request timed out (limit {}s)'.format(timeout)
        res.revoke()
        return JsonResponse(resp, status=408)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return JsonResponse(resp, status=400)

    resp['precursors'] = precursors
    for precursor in precursors:
        precursor['templates'] = precursor.pop('tforms')
    return JsonResponse(resp)
