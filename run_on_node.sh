#!/bin/bash

# specify resources
#SBATCH -G 1
#SBATCH -N 1
#SBATCH -e slurm-%j.log
#SBATCH -o slurm-%j.log

# max wallclock time
#SBATCH -t 48:00:00
# jobname
#SBATCH -J Train_wo_maps

# queue
#SBATCH -p lindahl5 -w erco-gpu01

module load cuda/12
source setup_chimerax.sh
# python3 data_generation/generate_input_maps.py -s 16400 -e 17807 -p 18
# python3 data_generation/generate_target_maps.py -s 17000 -e 17807 -p 22
# python3 data_generation/generate_dataset.py -p 22
# python3 data_generation/compute_rmsd_and_bonds.py -p 22
python3 model/generate_good_maps_with_model.py
# python3 testing/test_for_experimental_maps.py