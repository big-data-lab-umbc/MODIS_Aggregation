#!/bin/bash
#SBATCH --job-name=MODIS_agg
#SBATCH --output=modis_files_v5.out
#SBATCH --error=modis_files_v5.err
#SBATCH --partition=high_mem
#SBATCH --qos=medium+
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=MaxMemPerNode
#SBATCH --account=pi_zzbatmos

# Guide of running the script:
# python aggre_stats_series_v5.py <1/0 for Minimum> <1/0 for Maximum> <1/0 for Mean> <1/0 for Pixel Counts>
#                                 <1/0 for Standard Deviation> <1/0 for Histogram Counts> <1/0 for Joint Histogram>
#                                 <Variable Name List with Histogram Interval> 
#								  <Variable Name List with Joint Histogram Interval>  

# Example 1: Calculate Minimum, Maximum, Mean, Pixel Counts, Standard Deviation, Histogram and Joint Histogram. 
#			 Variables in input_file.csv & input_Jhist.csv:  
#              cloud_top_pressure_1km   
#              cloud_top_temperature_1km
#              cloud_emissivity_1km     
#              Cloud_Fraction           
#              cloud_top_height_1km     


#python aggre_stats_series_v2.py 1 1 1 1 1 0 1 filename.dat 21 150 320 histname.dat 4 0 1100
python3 aggre_stats_series_v5.py 1 1 1 1 1 1 1 input_file.csv input_Jhist.csv

