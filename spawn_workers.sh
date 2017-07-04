source activate askcos
export CELERY_LOG_LEVEL="WARN"

# Purge reservable queues to avoid weird conflicts
celery amqp queue.purge tb_worker_reservable
celery amqp queue.purge cr_worker_reservable
celery amqp queue.purge tb_worker
celery amqp queue.purge fp_worker
celery amqp queue.purge cr_worker

# The core workers 
# note: tb_worker is avilable for one-step even if tb_coordinator reserves all others
celery -A askcos_site worker -c 2 -Q tb_worker -n "treebuilder_worker@$(hostname)" --max-tasks-per-child 10000 --loglevel=$CELERY_LOG_LEVEL --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 2 -Q tb_coordinator -n "treebuilder_coordinator@$(hostname)"  --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 6 -Q fp_worker -n "forwardpredictor_worker@$(hostname)" --max-tasks-per-child 10000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 2 -Q fp_coordinator -n "forwardpredictor_coordinator@$(hostname)" --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 2 -Q context_worker -n "context_rec@$(hostname)" --max-tasks-per-child 50000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 6 -Q cr_worker -n "chiralretro_worker@$(hostname)" --max-tasks-per-child 10000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 2 -Q cr_coordinator -n "chiralretro_coordinator@$(hostname)" --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &

# Some bonus worker pools for the treebuilder(s) to reserve
celery -A askcos_site worker -c 6 -Q tb_worker,tb_worker_reservable -n "treebuilder_worker_pool1@$(hostname)" --max-tasks-per-child 10000 --loglevel=$CELERY_LOG_LEVEL --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 6 -Q tb_worker,tb_worker_reservable -n "treebuilder_worker_pool2@$(hostname)" --max-tasks-per-child 10000 --loglevel=$CELERY_LOG_LEVEL --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 6 -Q cr_worker,cr_worker_reservable -n "chiralretro_worker_pool1@$(hostname)" --max-tasks-per-child 10000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 6 -Q cr_worker,cr_worker_reservable -n "chiralretro_worker_pool2@$(hostname)" --max-tasks-per-child 10000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &


source deactivate