source activate askcos
celery -A askcos_site worker -c 4 -Q tb_worker -n "treebuilder_worker" &
celery -A askcos_site worker -c 1 -Q tb_coordinator -n "treebuilder_coordinator" &
celery -A askcos_site worker -c 4 -Q fp_worker -n "forwardpredictor_worker" &
celery -A askcos_site worker -c 1 -Q fp_coordinator -n "forwardpredictor_coordinator" &
celery -A askcos_site worker -c 1 -Q context_worker -n "context_rec" &
source deactivate 