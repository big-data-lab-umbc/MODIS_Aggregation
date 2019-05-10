#!/bin/bash

#SBATCH --job-name=cmacMODAgg
#SBATCH --output=slurm-cmacMODAgg_1Day_2vars.out
#SBATCH --partition=batch
#SBATCH --qos=medium+
#SBATCH --mem=32000
python Zhibo_MODIS_daily_mean.py MODAgg_out_2vars_20080101.hdf5 
