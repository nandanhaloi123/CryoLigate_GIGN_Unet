#!/bin/bash

#SBATCH --gpus 4

#SBATCH -t 3-00:00:00

#SBATCH -J ModelTrain

#SBATCH -e slurm-%j.log

#SBATCH -o slurm-%j.log

python Training/train_Unet3D_transformer_VAE.py

EOF
