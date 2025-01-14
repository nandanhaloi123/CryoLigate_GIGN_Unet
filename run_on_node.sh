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
#SBATCH -p lindahl5 -w erco-gpu02

module load cuda/12
source setup_chimerax.sh
# python3 cryoEM_maps/main.py -s 16400 -e 17807 -p 18
# python3 cryoEM_maps/main.py -s 16000 -e 17807 -p 18
# python3 generate_text_files_for_dataset.py -p 22
# python3 good_resolution_maps.py -s 17000 -e 17807 -p 22
# python3 preprocessing.py -s 0 -e 7565 -p 22
# python3 dataset_GIGN.py -p 22
# python3 train_example.py
python3 compute_rmsd.py -p 22
# python3 train_1.py
# python3 generate_good_maps_with_model.py
# python3 test_for_experimental_maps.py