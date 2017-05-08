source activate askcos
celery -A askcos_site worker -c 4 -Q tb_worker -l info -n "treebuilder_worker" &
celery -A askcos_site worker -c 1 -Q tb_coordinator -l info -n "treebuilder_coordinator" &
celery -A askcos_site worker -c 4 -Q fp_worker -l info -n "forwardpredictor_worker" &
celery -A askcos_site worker -c 1 -Q fp_coordinator -l info -n "forwardpredictor_coordinator" &
celery -A askcos_site worker -c 1 -Q context_worker -l info -n "context_rec" &
source deactivate 