#!/bin/bash

#SBATCH --job-name=cmacMODAgg
#SBATCH --output=slurm-cmacMODAgg.out
#SBATCH --partition=batch
#SBATCH --qos=medium+
#SBATCH --mem=32000
python Zhibo_MODIS_daily_mean.py MODAgg_out.hdf5 
