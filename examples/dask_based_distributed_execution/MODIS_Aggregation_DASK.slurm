#!/bin/bash
#SBATCH --partition=high_mem
#SBATCH --exclusive
#SBATCH --qos=medium+
#SBATCH --output=slurm-%x-%j-%u.out
#SBATCH --error=slurm-%x-%j-%u.out
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --account=pi_jianwu

# Example 1: Choose the time from 2008/01/01 to 2008/01/01 (Daily) with the region [-90,90,-180,180] with grid size [1,1] with sampling rate 1
#			 Calculate Minimum, Maximum, Mean, Pixel Counts, Standard Deviation, Histogram and Joint Histogram.
#			 Variables (5km resolution) in input_file_1km.csv & input_Jhist_1km.csv:
#              Cloud_Top_Pressure        
#              Cloud_Top_Temperature     
#              Cloud_Effective_Emissivity
#              cloud_fraction_CM         
#              Cloud_Top_Height          
python3 MODIS_Aggregation_DASK.py ../data_path.csv 2008/01/01 2008/01/01 [-90,90,-180,180] [1,1] [1] 1 1 1 1 1 1 1 ../input_file_5km.csv ../input_Jhist_5km.csv

# Example 2: Choose the time from 2008/01/01 to 2008/01/31 (Monthly) with the region [-90,90,-180,180] with grid size [1,1] with sampling rate 1
#			 Calculate Minimum, Maximum, Mean, Pixel Counts, Standard Deviation, Histogram and Joint Histogram.
#			 Variables (5km resolution) in input_file_5km.csv & input_Jhist_5km.csv:
#              Cloud_Top_Pressure        
#              Cloud_Top_Temperature     
#              Cloud_Effective_Emissivity
#              cloud_fraction_CM         
#              Cloud_Top_Height          
#python3 dask_based_modis_aggregation.py ../data_path.csv 2008/01/01 2008/01/31 [-90,90,-180,180] [1,1] [1] 1 1 1 1 1 1 1 ../input_file_5km.csv ../input_Jhist_5km.csv

# Example 3: Choose the time from 2008/01/01 to 2008/01/16 with the region [-20,20,0,35] with grid size [1,1] with sampling rate 2
#			 Calculate Minimum, Maximum, Mean, Pixel Counts.
#			 Variables (1km resolution) in input_file_1km.csv & input_Jhist_1km.csv:
#              cloud_top_pressure_1km   
#              cloud_top_temperature_1km
#              cloud_emissivity_1km     
#              cloud_fraction_CM        
#              cloud_top_height_1km     
#python3 dask_based_modis_aggregation.py ../data_path.csv 2008/01/01 2008/01/16 [-20,20,0,35] [1,1] [5] 1 1 1 0 0 0 0 ../input_file_1km.csv ../input_Jhist_1km.csv

# Example 4: Choose the time from 2008/01/01 to 2008/01/09 with the region [-90,90,-180,180] with grid size [0.5,0.625] with sampling rate 3
#			 Calculate Standard Deviation, Histogram and Joint Histogram.
#			 Variables (1km resolution) in input_file_1km.csv & input_Jhist_1km.csv:
#              cloud_top_pressure_1km   
#              cloud_top_temperature_1km
#              cloud_emissivity_1km     
#              cloud_fraction_CM        
#              cloud_top_height_1km     
#python3 dask_based_modis_aggregation.py ../data_path.csv 2008/01/01 2008/01/09 [-90,90,-180,180] [0.5,0.625] [5] 0 0 0 0 1 1 1 ../input_file_1km.csv ../input_Jhist_1km.csv
