#!/bin/bash
#SBATCH --job-name=mpi_Job
#SBATCH --partition=high_mem
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --qos=normal+
#SBATCH --mem=320GB
#SBATCH --output=slurm-%x-%j-%u.out
#SBATCH --error=slurm-%x-%j-%u.out

mpirun -N 90 python3 monthly-aggregation-file-level-mpi.py
