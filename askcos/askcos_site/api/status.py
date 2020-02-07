from django.http import JsonResponse
from askcos_site.celery import app

READABLE_NAMES = {
    'cr_network_worker': 'Context Recommender Worker',
    'ft_worker': 'Forward Predictor Worker',
    'sc_coordinator': 'Forward Predictor Scoring Coordinator',
    'tb_c_worker': 'One-Step/Tree Builder Retrosynthesis Worker',
    'tb_c_worker_preload': 'One-Step/Tree Builder Retrosynthesis Worker (Pre-loaded)',
    'tb_coordinator_mcts': 'Tree Builder Coordinator',
    'te_coordinator': 'Tree Evaluation Coordinator',
    'sites_worker': 'Site Selectivity Worker',
    'impurity_worker': 'Impurity worker',
    'atom_mapping_worker': 'Atom mapping worker',
}

def celery_status(request):
    resp = {}
    status = {}
    stats = app.control.inspect().stats()
    active = app.control.inspect().active()
    if not stats or not active:
        return JsonResponse(resp)
    worker_names = stats.keys()
    for worker in worker_names:
        name, server = worker.split('@')
        if not status.get(name):
            status[name] = {'available': 0, 'busy': 0}
        status[name]['busy'] += len(active[worker])
        status[name]['available'] += stats[worker]['pool']['max-concurrency'] - status[name]['busy']
    status_list = []
    for key in status:
        status_list.append({
            'name': READABLE_NAMES.get(key),
            'queue': key,
            'busy': status[key]['busy'],
            'available': status[key]['available']
        })
    for key, val in READABLE_NAMES.items():
        if key not in status:
            status_list.append({
                'name': READABLE_NAMES.get(key),
                'queue': key,
                'busy': 0,
                'available': 0
            })
    resp['queues'] = sorted(status_list, key=lambda x: x['name'])
    return JsonResponse(resp)
