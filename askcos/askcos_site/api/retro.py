from rdkit import Chem
from django.http import JsonResponse
from celery.exceptions import TimeoutError
from makeit import global_config as gc
from askcos_site.askcos_celery.treebuilder.tb_c_worker import get_top_precursors as get_top_precursors_c
from askcos_site.askcos_celery.treebuilder.tb_c_worker_preload import get_top_precursors as get_top_precursors_p

TIMEOUT = 120

def singlestep(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    run_async = request.GET.get('async', False)
    target = request.GET.get('target')
    max_num_templates = int(request.GET.get('num_templates', 100))
    max_cum_prob = float(request.GET.get('max_cum_prob', 0.995))
    fast_filter_threshold = float(request.GET.get('filter_threshold', 0.75))
    template_set = request.GET.get('template_set', 'reaxys')
    template_prioritizer = request.GET.get('template_prioritizer', 'reaxys')

    if template_set not in ['reaxys', 'uspto']:
        resp['error'] = 'Template set {} not available'.format(template_set)
        return JsonResponse(resp, status=400)

    if template_prioritizer not in ['reaxys', 'uspto', 'uspto_pretrained']:
        resp['error'] = 'Template prioritizer {} not available'.format(template_prioritizer)
        return JsonResponse(resp, status=400)

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

    if max_cum_prob > 0.999 and max_num_templates > 1000:
        res = get_top_precursors_p.delay(
            target,
            fast_filter_threshold=fast_filter_threshold,
            max_cum_prob=max_cum_prob,
            max_num_templates=max_num_templates,
            cluster=cluster,
            cluster_method=cluster_method,
            cluster_feature=cluster_feature,
            cluster_fp_type=cluster_fp_type,
            cluster_fp_length=cluster_fp_length,
            cluster_fp_radius=cluster_fp_radius
        )
    else:
        res = get_top_precursors_c.delay(
            target,
            template_set=template_set,
            template_prioritizer=template_prioritizer,
            fast_filter_threshold=fast_filter_threshold,
            max_cum_prob=max_cum_prob,
            max_num_templates=max_num_templates,
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

    try:
        (smiles, precursors) = res.get(TIMEOUT)
    except TimeoutError:
        resp['error'] = 'API request timed out (limit {}s)'.format(TIMEOUT)
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
