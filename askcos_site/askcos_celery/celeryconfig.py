TASK_SERIALIZER = 'json'
RESULT_SERIALIZER = 'json'
ACCEPT_CONTENT = ['json']
TIMEZONE = 'US/Eastern'
ENABLE_UTC = True

RESULT_EXPIRES = 600 # only keep results for 10 minutes max

TASK_ROUTES = {
    'askcos_site.askcos_celery.treebuilder.worker.*': {'queue': 'tb_worker'},
    'askcos_site.askcos_celery.treebuilder.coordinator.*': {'queue': 'tb_coordinator'},
    'askcos_site.askcos_celery.forwardpredictor.worker.*': {'queue': 'fp_worker'},
    'askcos_site.askcos_celery.forwardpredictor.coordinator.*': {'queue': 'fp_coordinator'}, 
    'askcos_site.askcos_celery.contextrecommender.worker.*': {'queue': 'context_worker'},   
}

CELERY_ROUTES = TASK_ROUTES