#!/bin/bash
export PATH="/home/mpiuser/miniconda2/bin:$PATH"
export PYTHONPATH="/home/mpiuser/miniconda2/envs/askcos/lib/python2.7/site-packages:/home/mpiuser/ASKCOS/Make-It:/home/mpiuser/ASKCOS/ASKCOS_Website:$PYTHONPATH"
echo $PATH

source activate askcos
export CELERY_LOG_LEVEL="WARN"

# Tree evaluator coordinator
celery -A askcos_site worker -c 2 -Q te_coordinator -n "te_coordinator@$(hostname)"  --loglevel=$CELERY_LOG_LEVEL --maxtasksperchild 50 --logfile=/home/mpiuser/ASKCOS/ASKCOS_Website/celery_logs/%p.log &

# Scoring coordinator
celery -A askcos_site worker -c 4 -Q sc_coordinator -n "sc_coordinator@$(hostname)"  --loglevel=$CELERY_LOG_LEVEL --maxtasksperchild 1000 --logfile=/home/mpiuser/ASKCOS/ASKCOS_Website/celery_logs/%p.log &

# Template-based forward predictor
celery -A askcos_site worker -c 6 -Q ft_worker -n "ft_worker@$(hostname)" --loglevel=$CELERY_LOG_LEVEL  --max-tasks-per-child 10000 --loglevel=$CELERY_LOG_LEVEL  --logfile=/home/mpiuser/ASKCOS/ASKCOS_Website/celery_logs/%p.log &

# Template-free forward predictor (for testing)
# celery -A askcos_site worker -c 1 -Q tffp_worker -n "tffp_worker@$(hostname)" --max-tasks-per-child 10000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &

# Context recommender coordinator
celery -A askcos_site worker -c 4 -Q cr_coordinator -n "cr_coordinator@$(hostname)" --loglevel=$CELERY_LOG_LEVEL --max-tasks-per-child 10000 --loglevel=$CELERY_LOG_LEVEL  --logfile=/home/mpiuser/ASKCOS/ASKCOS_Website/celery_logs/%p.log &

# Context recommender nearest-neighbor worker
# celery -A askcos_site worker -c 2 -Q cr_nn_worker -n "cr_nn_worker@$(hostname)" --max-tasks-per-child 1000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &

# Context recommender neural network worker
celery -A askcos_site worker -c 4 -Q cr_network_worker -n "cr_network_worker@$(hostname)" --loglevel=$CELERY_LOG_LEVEL --max-tasks-per-child 50000 --loglevel=$CELERY_LOG_LEVEL  --logfile=/home/mpiuser/ASKCOS/ASKCOS_Website/celery_logs/%p.log &

source deactivate