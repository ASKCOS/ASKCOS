source activate askcos
celery -A askcos_site worker -c 4 -Q tb_worker -n "treebuilder_worker@$(hostname)" --max-tasks-per-child 10000 --loglevel=INFO --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 1 -Q tb_coordinator -n "treebuilder_coordinator@$(hostname)"  --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 4 -Q fp_worker -n "forwardpredictor_worker@$(hostname)" --max-tasks-per-child 10000 --loglevel=INFO  --logfile=%p.log &
celery -A askcos_site worker -c 1 -Q fp_coordinator -n "forwardpredictor_coordinator@$(hostname)" --loglevel=INFO  --logfile=%p.log &
celery -A askcos_site worker -c 1 -Q context_worker -n "context_rec@$(hostname)" --max-tasks-per-child 50000 --loglevel=INFO  --logfile=%p.log &
source deactivate