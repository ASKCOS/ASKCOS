source activate askcos
export CELERY_LOG_LEVEL="WARN"
rm celery_logs/*.log

# Purge reservable queues to avoid weird conflicts
celery amqp queue.purge tb_worker_reservable
celery amqp queue.purge cr_worker_reservable
celery amqp queue.purge tb_worker
celery amqp queue.purge fp_worker
celery amqp queue.purge cr_worker


# The core workers 
# note: tb_worker is avilable for one-step even if tb_coordinator reserves all others
celery -A askcos_site worker -c 4 -Q ft_worker -n "forward_transformer_worker_pool1@$(hostname)" --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 4 -Q tb_worker,tb_worker_reservable -n "treebuilder_worker_pool1@$(hostname)" --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 4 -Q tb_c_worker, tb_c_worker_reservable -n "treebuilder_chiral_worker_pool1@$(hostname)" --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 1 -Q tb_coordinator -n "treebuilder_coordinator_pool1@$(hostname)" --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 1 -Q te_coordinator -n "treeevaluator_coordinator_pool1@$(hostname)" --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 1 -Q sc_coordinator -n "scoring_coordinator_pool1@$(hostname)" --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 1 -Q cr_nn_worker -n "neareast_neighbor_context_rec_pool1@$(hostname)" --logfile=celery_logs/%p.log &
celery -A askcos_site worker -c 1 -Q cr_coordinator -n "context_recommender_coordinator_pool1@$(hostname)" --logfile=celery_logs/%p.log &


source deactivate