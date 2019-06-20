#!/bin/bash

#SBATCH --job-name=cmacMODAgg
#SBATCH --output=slurm-MODAgg_3var_parMonth_200801.out
#SBATCH --partition=batch
#SBATCH --qos=medium+
#SBATCH --mem=8000
#SBATCH --array=1-31%31
mode='parMonth'
day=$SLURM_ARRAY_TASK_ID
python MODAgg_daily_mean.py MODAgg_3var_parMonth_ ${mode} ${day}
