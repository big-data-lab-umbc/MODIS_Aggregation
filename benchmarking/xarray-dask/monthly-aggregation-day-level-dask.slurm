#!/bin/bash -l
#SBATCH --job-name=mod_agg
#SBATCH --output=%x-%j_dask-monthly-aggregation-day-level.out
#SBATCH --error=%x-%j_dask-monthly-aggregation-day-level.err
#SBATCH --partition=batch
#SBATCH --qos=long+
#SBATCH --mem=MaxMemPerNode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=pi_jianwu
#SBATCH --exclusive

python3 monthly-aggregation-day-level-dask.py

