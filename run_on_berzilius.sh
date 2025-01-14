#!/bin/bash

#SBATCH --gpus 1

#SBATCH -t 3-00:00:00

#SBATCH -J ModelTrain

#SBATCH -e slurm-%j.log

#SBATCH -o slurm-%j.log

python train_example_berzilius.py

EOF