#!/bin/bash
#SBATCH --job-name=mpi_Job-file-level
#SBATCH --partition=high_mem
##SBATCH --nodes=5
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --qos=normal+
#SBATCH --mem=320GB
#SBATCH --output=slurm-%x-%j-%u-%N.out
#SBATCH --error=slurm-%x-%j-%u-%N.out

##sbatch command example: sbatch -N2 runMPI-monthly-aggregation-day-level.slurm

echo "node number: $SLURM_NNODES"
echo "nodes: $(echo $SLURM_NODELIST | cut -d'[' -f 2 | cut -d']' -f 1)"
mpirun -N $((16*$SLURM_NNODES)) python3 monthly-aggregation-file-level-mpi.py
