#!/bin/bash

#SBATCH --gpus 1

#SBATCH -t 3-00:00:00

#SBATCH -J SCUnet3D

#SBATCH -e slurm-%j.log

#SBATCH -o slurm-%j.log
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python Training/train_SCUnet3D_L1_Loss.py --alpha 20 --beta 0.1
# torchrun --nproc_per_node=2 Training/train_Unet3D_transformer_VAE_GAN.py

# torchrun --nproc_per_node=8 Training/train_Unet3D_transformer_VAE_GAN.py

