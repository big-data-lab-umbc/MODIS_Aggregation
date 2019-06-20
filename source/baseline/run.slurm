#!/bin/bash

#SBATCH --job-name=cmacMODAgg
#SBATCH --output=slurm-MODAgg_3var_month_200801.out
#SBATCH --partition=batch
#SBATCH --qos=long+
#SBATCH --mem=32000
mode='month'
day=1
python MODAgg_daily_mean.py MODAgg_3var_month_ ${mode} ${day}
