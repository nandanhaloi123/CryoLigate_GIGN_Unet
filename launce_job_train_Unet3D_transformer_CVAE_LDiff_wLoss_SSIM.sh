#!/bin/bash

# Create 20 combinations of alpha, beta, gamma
alphas=(1.0 0.5 0.1 2.0 3.0 4.0 5.0)
betas=(0.001 0.01 0.1 0.5 1.0)
gammas=(0.001 0.01 0.1 0.5 1.0)

job_count=0

for alpha in "${alphas[@]}"; do
for alpha2 in "${alphas[@]}"; do
  for beta in "${betas[@]}"; do
    for gamma in "${gammas[@]}"; do
      if [ $job_count -ge 50 ]; then
        break 3  # Exit all loops after 30 jobs
      fi

      job_script="submit_alpha_${alpha}_alpha2_${alpha2}_beta_${beta}_gamma_${gamma}.sh"

      cat <<EOF > $job_script
#!/bin/bash
#SBATCH --gpus 1
#SBATCH -t 3-00:00:00
#SBATCH -J CVAED_${job_count}
#SBATCH -e slurm-%j.log
#SBATCH -o slurm-%j.log
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python Training/train_Unet3D_transformer_CVAE_LDiff_wLoss_SSIM.py --alpha ${alpha} --alpha2 ${alpha2} --beta ${beta} --gamma ${gamma}
EOF
      chmod +x $job_script
      sbatch $job_script

      ((job_count++))
    done
  done
done
done