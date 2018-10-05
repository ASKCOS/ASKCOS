#!/bin/bash
export PATH="/home/mpiuser/miniconda2/bin:$PATH"
export PYTHONPATH="/home/mpiuser/miniconda2/envs/askcos/lib/python2.7/site-packages:/home/mpiuser/ASKCOS/Make-It:/home/mpiuser/ASKCOS/ASKCOS_Website:$PYTHONPATH"
echo $PATH

source activate askcos
export CELERY_LOG_LEVEL="WARN"

# Coordinators
KERAS_BACKEND=theano THEANO_FLAGS=device=cpu celery -A askcos_site worker -c 2 -Q tb_coordinator_mcts -n "tb_coordinator_mcts@$(hostname)"  --maxtasksperchild 500 --logfile=celery_logs/%p.log &

# Builders
KERAS_BACKEND=theano THEANO_FLAGS=device=cpu celery -A askcos_site worker -c 18 -Q tb_c_worker -n "tb_c_worker_pool1@$(hostname)" --maxtasksperchild 10000  --logfile=celery_logs/%p.log &

source deactivate