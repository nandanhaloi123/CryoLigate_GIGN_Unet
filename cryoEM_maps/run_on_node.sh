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
#SBATCH -p lindahl5 -w erco-gpu03

module load cuda/12
source setup_chimerax.sh
python3 cryoEM_maps/main.py -s 2000 -e 2999 -p 18
