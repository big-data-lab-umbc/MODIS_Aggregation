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

# Example 1: Calculate Minimum, Maximum, Mean, Pixel Counts, Standard Deviation, Joint Histogram. 
#			 Variables in filename.dat: cloud_top_pressure, cloud_top_temperature, Cloud_Optical_Thickness.
#			 Bin number for histogram counts is 21. 
#			 Boundary for histogram is 150 & 320 (cloud_top_temperature)
#			 Variables in histname.dat: cloud_top_temperature, cloud_top_pressure
#			 Bin number for the second variable for joint histrogram is 4.
#			 Boundary for the second variable for joint histogram is 0 & 1100 (cloud_top_pressure).

python aggre_stats_series_v2.py 1 1 1 1 1 0 1 filename.dat 21 150 320 histname.dat 4 0 1100 



