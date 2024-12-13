#!/bin/bash

# specify resources
#SBATCH -G 1
#SBATCH -N 1
#SBATCH -e slurm-%j.log
#SBATCH -o slurm-%j.log

# max wallclock time
#SBATCH -t 48:00:00
# jobname
#SBATCH -J Complexes

# queue
#SBATCH -p lindahl5 -w erco-gpu01

module load cuda/12
source setup_chimerax.sh
python3 preprocessing.py -s 0 -e 7565 -p 22
