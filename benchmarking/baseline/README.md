#Demo command to run the same slurm file with different sampling rates
sampling=3 && sbatch --export=sampling=$sampling --job-name=$sampling-sampling modis-baseline-job-file-level-sampling.slurm
