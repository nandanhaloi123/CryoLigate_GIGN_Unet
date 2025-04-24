#!/bin/bash

# Define your hyperparameter values
alphas=(1.0 3.0 5.0)
betas=(1.0 3.0 5.0)

# Optional: limit the number of jobs
MAX_JOBS=-1  # Set to -1 to run all combinations
job_count=0

# Loop over all combinations
for alpha in "${alphas[@]}"; do
    for beta in "${betas[@]}"; do
        # Check limit
        if [[ $MAX_JOBS -ge 0 && $job_count -ge $MAX_JOBS ]]; then
          echo "Reached max job limit: $MAX_JOBS"
          exit 0
        fi

        # Create job script name
        job_script="submit_alpha_${alpha}_beta_${beta}.sh"

        # Write job script
        cat <<EOF > $job_script
#!/bin/bash
#SBATCH --gpus 1
#SBATCH -t 3-00:00:00
#SBATCH -J Unet${job_count}
#SBATCH -e slurm-%j.log
#SBATCH -o slurm-%j.log
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python Training/train_Unet3D_Enc_Embedding.py --alpha ${alpha} --beta ${beta}
EOF

        chmod +x $job_script
        sbatch $job_script
        ((job_count++))
        rm $job_script 
      done
    done

echo "Submitted $job_count jobs."
