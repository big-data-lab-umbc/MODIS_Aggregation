#!/bin/bash
#SBATCH --partition=high_mem
#SBATCH --exclusive
#SBATCH --qos=medium+
#SBATCH --output=slurm-%x-%j-%u.out
#SBATCH --error=slurm-%x-%j-%u.out
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --account=pi_jianwu

python3 /umbc/xfs1/jianwu/common/MODIS_Aggregation/xin-code/MODIS_Aggregation/source/baseline/MODIS_Aggregation_dask_chunk_jianwu_local_filepair_xin_n8.py data_path.csv 2008/01/01 2008/01/31 [-90,90,-180,180] [1,1] [5] 1 1 1 1 1 1 1 input_file_1km.csv input_Jhist_1km.csv
