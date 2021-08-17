#!/bin/bash
#SBATCH --job-name=MODISMPI
#SBATCH --output=journal_benchmark/MPI/modis_MPI_5km.out
#SBATCH --error=journal_benchmark/MPI/modis_MPI_5km.err
#SBATCH --partition=high_mem
#SBATCH --qos=medium+
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=20
##SBATCH --exclusive
#SBATCH --time=12:00:00
#SBATCH --mem=MaxMemPerNode
#SBATCH --account=pi_jianwu

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
#			 Variables (1km resolution) in input_file_1km.csv & input_Jhist_1km.csv:
#              cloud_top_pressure_1km
#              cloud_top_temperature_1km
#              cloud_emissivity_1km
#              Cloud_Fraction
#              cloud_top_height_1km
mpirun -n 80 MODIS_Aggregation_MPI.py ../data_path.csv 2008/01/01 2008/01/31 [-90,90,-180,180] [1,1] [5] 1 1 1 1 1 1 1 ../input_file_1km.csv ../input_Jhist_1km.csv

# Example 2: Choose the time from 2008/01/01 to 2008/01/31 with the region [-90,90,-180,180] with grid size [1,1] with sampling rate 3
#			 Calculate Minimum, Maximum, Mean, Pixel Counts, Standard Deviation, Histogram and Joint Histogram.
#			 Variables (5km resolution) in input_file_5km.csv & input_Jhist_5km.csv:
#              Cloud_Top_Pressure
#              Cloud_Top_Temperature
#              Cloud_Effective_Emissivity
#              cloud_fraction
#              Cloud_Top_Height
#mpirun -n 80 python3 MODIS_Aggregation_MPI.py data_path.csv 2008/01/01 2008/01/31 [-90,90,-180,180] [1,1] [5] 1 1 1 1 1 1 1 input_file_5km.csv input_Jhist_5km.csv

# Example 3: Choose the time from 2008/01/01 to 2008/01/16 with the region [-20,20,0,35] with grid size [1,1] with sampling rate 2
#			 Calculate Minimum, Maximum, Mean, Pixel Counts.
#			 Variables (1km resolution) in input_file_1km.csv & input_Jhist_1km.csv:
#              cloud_top_pressure_1km
#              cloud_top_temperature_1km
#              cloud_emissivity_1km
#              Cloud_Fraction
#              cloud_top_height_1km
#mpirun -n 80 python3 MODIS_Aggregation_MPI.py data_path.csv 2008/01/01 2008/01/16 [-20,20,0,35] [1,1] [5] 1 1 1 0 0 0 0 input_file_1km.csv input_Jhist_1km.csv

# Example 4: Choose the time from 2008/01/01 to 2008/01/09 with the region [-90,90,-180,180] with grid size [0.5,0.625] with sampling rate 3
#			 Calculate Standard Deviation, Histogram and Joint Histogram.
#			 Variables (1km resolution) in input_file_1km.csv & input_Jhist_1km.csv:
#              cloud_top_pressure_1km
#              cloud_top_temperature_1km
#              cloud_emissivity_1km
#              Cloud_Fraction
#              cloud_top_height_1km
#mpirun -n 80 python3 MODIS_Aggregation_MPI.py data_path.csv 2008/01/01 2008/01/09 [-90,90,-180,180] [0.5,0.625] [5] 0 0 0 0 1 1 1 input_file_1km.csv input_Jhist_1km.csv
