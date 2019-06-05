#!/bin/bash

#SBATCH --job-name=cmacMODAgg
#SBATCH --output=slurm-MODAgg_3var_oneDay_2008010.out
#SBATCH --partition=batch
#SBATCH --qos=medium+
#SBATCH --mem=32000
python Zhibo_MODIS_daily_mean.py MODAgg_3var_oneDay_20080101.hdf5 oneDay
