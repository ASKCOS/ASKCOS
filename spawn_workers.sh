source activate askcos
celery -A askcos_site worker -c 4 -Q tb_worker -n "treebuilder_worker" --max-tasks-per-child 10000 &
celery -A askcos_site worker -c 1 -Q tb_coordinator -n "treebuilder_coordinator" &
celery -A askcos_site worker -c 4 -Q fp_worker -n "forwardpredictor_worker" --max-tasks-per-child 10000 &
celery -A askcos_site worker -c 1 -Q fp_coordinator -n "forwardpredictor_coordinator" &
celery -A askcos_site worker -c 1 -Q context_worker -n "context_rec" --max-tasks-per-child 50000 &
source deactivate



# source activate askcos
# celery -A askcos_site multi start tb_coordinator@%h fp_coordinator@%h context_worker@%h tb_worker@%h fp_worker@%h -c:1-3 1 -c:4-5 4 -Q:1 tb_coordinatotr -Q:2 fp_coordinator -Q:3 context_worker -Q:4 tb_worker -Q:5 fp_worker
# source deactivate