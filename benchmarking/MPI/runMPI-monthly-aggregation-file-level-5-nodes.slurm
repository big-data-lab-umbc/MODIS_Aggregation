#!/bin/bash
#SBATCH --job-name=sparkJob
#SBATCH --partition=high_mem
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=18
#SBATCH --exclusive
#SBATCH --qos=normal+
#SBATCH --mem=20GB
#SBATCH --output=slurm-%x-%j-%u.out
#SBATCH --error=slurm-%x-%j-%u.out

mpirun -N 90 python3 monthly-aggregation-file-level-mpi.py
