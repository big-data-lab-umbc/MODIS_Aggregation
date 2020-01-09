#!/bin/bash
#SBATCH --job-name=MODIS_agg
#SBATCH --output=modis_files.out
#SBATCH --error=modis_files.err
#SBATCH --partition=batch
#SBATCH --qos=medium+
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=MaxMemPerNode
#SBATCH --account=pi_zzbatmos

# Guide of running the script:
# python aggre_stats_series_v2.py <1/0 for Minimum> <1/0 for Maximum> <1/0 for Mean> <1/0 for Pixel Counts>
#                                 <1/0 for Standard Deviation> <1/0 for Histogram Counts> <1/0 for Joint Histogram>
#                                 <Variable Name List> <Bin Number for Histogram Counts> 
#                                 <lower Boundary> <upper Boundary> 
#								  <Joint Histogram Name List> <Bin Number for Joint Histogram> 
#								  <lower Boundary for joint histogram> <upper Boundary for joint histogram> 

python aggre_stats_series_v2.py 1 1 1 1 1 0 1 filename.dat 21 150 320 histname.dat 4 0 1100

