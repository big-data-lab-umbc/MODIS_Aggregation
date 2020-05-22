# Demo command to run the same slurm file with different sampling rates on different nodes

running the same slurm on 2 nodes with sampling rate as 3: node=2 && sampling=3 && sbatch --export=sampling=$sampling,node=$node --job-name=$node-node-$sampling-sampling --nodes=$node runSparkAppOnTaki-monthly-aggregation-file-level-sampling.slurm

sampling=1 for no sampling
