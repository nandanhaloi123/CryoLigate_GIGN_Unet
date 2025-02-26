#!/bin/bash

#SBATCH --gpus 1

#SBATCH -t 3-00:00:00

#SBATCH -J VAE_GAN

#SBATCH -e slurm-%j.log

#SBATCH -o slurm-%j.log

python Training/train_Unet3D_transformer_VAE_GAN.py
# torchrun --nproc_per_node=2 Training/train_Unet3D_transformer_VAE_GAN.py

# torchrun --nproc_per_node=8 Training/train_Unet3D_transformer_VAE_GAN.py

