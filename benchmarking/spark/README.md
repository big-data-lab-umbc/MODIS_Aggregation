# Demo command to run the same slurm file with different sampling rates on different nodes

running the same slurm on 2 nodes with sampling rate as 4: node=2 && sampling=4 && sbatch --export=sampling=$sampling,node=$node --job-name=file-level-$node-node-$sampling-sampling --nodes=$node runSparkAppOnTaki-monthly-aggregation-file-level-sampling.slurm

running the same slurm on 2 nodes with sampling rate as 3: node=2 && sampling=4 && sbatch --export=sampling=$sampling,node=$node --job-name=pixel-level-$node-node-$sampling-sampling --nodes=$node runSparkAppOnTaki-monthly-aggregation-pixel-level-sampling.slurm

running the same slurm on 3 nodes with the default sampling rate (3): node=2 && sbatch --export=node=$node --job-name=day-level-$node-node-3-sampling --nodes=$node runSparkAppOnTaki-monthly-aggregation-day-level-sampling.slurm

Note: 

1. sampling=1 for no sampling
2. default sampling rate is 3 is sampling parameter is not sepecified
