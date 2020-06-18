#!/bin/bash
#SBATCH --job-name=MODIS_agg
#SBATCH --output=modis_files_v8_dask2.out
#SBATCH --error=modis_files_v8_dask2.err
#SBATCH --partition=high_mem
#SBATCH --qos=medium+
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=MaxMemPerNode
#SBATCH --account=pi_zzbatmos

# Guide of running the script:
# python aggre_stats_series_v5.py <Data Path file> <Start Date> <End Date>
#								  <Polygon boundaries: [lower lat, upper lat, left lon, right lon> 
# 								  <Lat & Lon Grid Size: [lat grid size, lon grid size]> 
#								  <Sampling number larger than 0> 
#								  <1/0 for Minimum> <1/0 for Maximum> <1/0 for Mean> <1/0 for Pixel Counts>
#								  <1/0 for Standard Deviation> <1/0 for Histogram Counts> <1/0 for Joint Histogram>
#								  <Variable Name List with Histogram Interval> 
#								  <Variable Name List with Joint Histogram Interval>  

# Example 1: Choose the time from 2008/01/01 to 2008/01/31 with the region [-90,90,-180,180] with grid size [1,1] with sampling rate 3
#			 Calculate Minimum, Maximum, Mean, Pixel Counts, Standard Deviation, Histogram and Joint Histogram. 
#			 Variables in input_file.csv & input_Jhist.csv:  
#              cloud_top_pressure_1km   
#              cloud_top_temperature_1km
#              cloud_emissivity_1km     
#              Cloud_Fraction           
#              cloud_top_height_1km     
python3 baseline_series_v8_dask2.py data_path.csv 2008/01/01 2008/01/31 [-90,90,-180,180] [1,1] [5] 1 1 1 1 1 1 1 input_file.csv input_Jhist.csv

# Example 2: Choose the time from 2008/01/01 to 2008/01/16 with the region [-20,20,0,35] with grid size [1,1] with sampling rate 2
#			 Calculate Minimum, Maximum, Mean, Pixel Counts. 
#			 Variables in input_file.csv & input_Jhist.csv:  
#              cloud_top_pressure_1km   
#              cloud_top_temperature_1km
#              cloud_emissivity_1km     
#              Cloud_Fraction           
#              cloud_top_height_1km     
#python3 baseline_series_v8_dask.py data_path.csv 2008/01/01 2008/01/01 [-20,20,0,35] [1,1] [5] 1 1 1 0 0 0 0 input_file.csv input_Jhist.csv

# Example 3: Choose the time from 2008/01/01 to 2008/01/09 with the region [-90,90,-180,180] with grid size [0.5,0.625] with sampling rate 3
#			 Calculate Standard Deviation, Histogram and Joint Histogram. 
#			 Variables in input_file.csv & input_Jhist.csv:  
#              cloud_top_pressure_1km   
#              cloud_top_temperature_1km
#              cloud_emissivity_1km     
#              Cloud_Fraction           
#              cloud_top_height_1km     
#python3 baseline_series_v8.py data_path.csv 2008/01/01 2008/01/10 [-90,90,-180,180] [0.5,0.625] [5] 0 0 0 0 1 1 1 input_file.csv input_Jhist.csv
