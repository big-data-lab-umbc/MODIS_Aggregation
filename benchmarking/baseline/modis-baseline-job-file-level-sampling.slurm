#!/bin/bash
#SBATCH --partition=high_mem
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --qos=long+
#SBATCH --output=slurm-%x-%j-%u.out
#SBATCH --error=slurm-%x-%j-%u.out

#Demo command to run the same slurm file with different sampling rates: sampling=3 && sbatch --export=sampling=$sampling --job-name=$sampling-sampling modis-baseline-job-file-level-sampling.slurm

export CONDA_PREFIX=/umbc/xfs1/cybertrn/common/Softwares/anaconda3/

$CONDA_PREFIX/bin/python /umbc/xfs1/jianwu/common/MODIS_Aggregation/jianwu-code/github/MODIS_Aggregation/benchmarking/baseline/monthly-aggregation-file-level-for-loop-netcdf.py $sampling
