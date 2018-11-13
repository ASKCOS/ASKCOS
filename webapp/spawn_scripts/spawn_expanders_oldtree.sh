#!/bin/bash
export PATH="/home/mpiuser/miniconda2/bin:$PATH"
export PYTHONPATH="/home/mpiuser/miniconda2/envs/askcos/lib/python2.7/site-packages:/home/mpiuser/ASKCOS/Make-It:/home/mpiuser/ASKCOS/ASKCOS_Website:$PYTHONPATH"
echo $PATH

source activate askcos
export CELERY_LOG_LEVEL="WARN"

# Coordinators
KERAS_BACKEND=theano THEANO_FLAGS=device=cpu celery -A askcos_site worker -c 4 -Q tb_coordinator -n "tb_coordinator@$(hostname)"  --maxtasksperchild 500 --logfile=celery_logs/%p.log &

# Builders
KERAS_BACKEND=theano THEANO_FLAGS=device=cpu celery -A askcos_site worker -c 8 -Q tb_c_worker,tb_c_worker_reservable -n "tb_c_worker_pool_reservable1@$(hostname)" --maxtasksperchild 5000   --logfile=celery_logs/%p.log &
KERAS_BACKEND=theano THEANO_FLAGS=device=cpu celery -A askcos_site worker -c 8 -Q tb_c_worker,tb_c_worker_reservable -n "tb_c_worker_pool_reservable2@$(hostname)" --maxtasksperchild 5000  --logfile=celery_logs/%p.log &

source deactivate