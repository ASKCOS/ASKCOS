#!/bin/bash
#SBATCH --job-name=expansion_worker
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3GB

source activate askcos
export CELERY_LOG_LEVEL="WARN"

# Chiral retro transformers
celery -A askcos_site worker -c 8 -Q tb_c_worker,tb_c_worker_reservable -n "tb_c_worker_pool1@$(hostname)" --maxtasksperchild 5000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p_$RANDOM.log

source deactivate
