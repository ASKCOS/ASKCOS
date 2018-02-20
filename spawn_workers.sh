source activate askcos
export CELERY_LOG_LEVEL="WARN"
rm celery_logs/*.log

# Purge reservable queues to avoid weird conflicts
celery amqp queue.purge tb_worker
celery amqp queue.purge tb_worker_reservable
celery amqp queue.purge tb_c_worker
celery amqp queue.purge tb_c_worker_reservable
celery amqp queue.purge tb_coordinator
celery amqp queue.purge te_coordinator
celery amqp queue.purge sc_coordinator
celery amqp queue.purge ft_worker
celery amqp queue.purge tffp_worker
celery amqp queue.purge cr_coordinator
celery amqp queue.purge cr_nn_worker


# Achiral retro transformers
celery -A askcos_site worker -c 6 -Q tb_worker,tb_worker_reservable -n "tb_worker_pool1@$(hostname)" --maxtasksperchild 5000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &

# Chiral retro transformers
celery -A askcos_site worker -c 8 -Q tb_c_worker,tb_c_worker_reservable -n "tb_c_worker_pool1@$(hostname)" --maxtasksperchild 5000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 8 -Q tb_c_worker,tb_c_worker_reservable -n "tb_c_worker_pool2@$(hostname)" --maxtasksperchild 5000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 8 -Q tb_c_worker,tb_c_worker_reservable -n "tb_c_worker_pool3@$(hostname)" --maxtasksperchild 5000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 3 -Q tb_c_worker -n "tb_c_worker@$(hostname)" --maxtasksperchild 1000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &

# Tree builder coordinator
celery -A askcos_site worker -c 3 -Q tb_coordinator -n "tb_coordinator@$(hostname)"  --maxtasksperchild 5 --logfile=celery_logs/%p.log &

# Tree evaluator coordinator
celery -A askcos_site worker -c 2 -Q te_coordinator -n "te_coordinator@$(hostname)"  --maxtasksperchild 5 --logfile=celery_logs/%p.log &

# Scoring coordinator
celery -A askcos_site worker -c 2 -Q sc_coordinator -n "sc_coordinator@$(hostname)"  --maxtasksperchild 100 --logfile=celery_logs/%p.log &

# Template-based forward predictor
celery -A askcos_site worker -c 6 -Q ft_worker -n "ft_worker@$(hostname)" --max-tasks-per-child 10000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &

# Template-free forward predictor (for testing)
# celery -A askcos_site worker -c 1 -Q tffp_worker -n "tffp_worker@$(hostname)" --max-tasks-per-child 10000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &

# Context recommender coordinator
celery -A askcos_site worker -c 2 -Q cr_coordinator -n "cr_coordinator@$(hostname)" --max-tasks-per-child 1000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &

# Context recommender nearest-neighbor worker
celery -A askcos_site worker -c 2 -Q cr_nn_worker -n "cr_nn_worker@$(hostname)" --max-tasks-per-child 1000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &
# Context recommender neural network worker
celery -A askcos_site worker -c 2 -Q cr_network_worker -n "cr_network_worker@$(hostname)" --max-tasks-per-child 5000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &


source deactivate
