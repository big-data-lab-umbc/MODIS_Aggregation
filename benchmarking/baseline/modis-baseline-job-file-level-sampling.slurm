#!/bin/bash
#SBATCH --partition=high_mem
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --qos=long+
#SBATCH --output=slurm-%x-%j-%u.out
#SBATCH --error=slurm-%x-%j-%u.out

export CONDA_PREFIX=/umbc/xfs1/cybertrn/common/Softwares/anaconda3/

$CONDA_PREFIX/bin/python /umbc/xfs1/jianwu/common/MODIS_Aggregation/jianwu-code/github/MODIS_Aggregation/benchmarking/baseline/monthly-aggregation-file-level-for-loop-netcdf.py $sampling
