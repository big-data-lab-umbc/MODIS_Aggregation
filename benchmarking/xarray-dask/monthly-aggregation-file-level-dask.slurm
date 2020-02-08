#!/bin/bash -l
#SBATCH --job-name=mod_agg
#SBATCH --output=%x-%j_dask-monthly-aggregation-file-level.out
#SBATCH --error=%x-%j_dask-monthly-aggregation-file-level.err
#SBATCH --partition=batch
#SBATCH --qos=long+
#SBATCH --mem=MaxMemPerNode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=pi_jianwu
#SBATCH --exclusive

srun /umbc/xfs1/cybertrn/common/Softwares/anaconda3/bin/python monthly-aggregation-file-level-dask.py
