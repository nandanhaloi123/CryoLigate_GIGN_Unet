#!/bin/bash

# specify resources
#SBATCH -G 1
#SBATCH -N 1
#SBATCH -e slurm-%j.log
#SBATCH -o slurm-%j.log

# max wallclock time
#SBATCH -t 48:00:00
# jobname
#SBATCH -J Maps

# queue
#SBATCH -p lindahl5 -w erco-gpu17

module load cuda/12
source setup_chimerax.sh
python3 good_resolution_maps.py -s 6000 -e 7565 -p 18
