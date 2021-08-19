#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Main Program: Run MODIS AGGREGATION IN SERIES WITH FLEXIBLE STATISTICS

Created on 2019

@author: Jianyu Zheng

V2 Updates: Add statistics for flexible variables

V3 Updates: Add 1d histogram with upper and lower boundaries

V4 Updates: Add 2d histogram by using additional input file

V5 Updates: Refine 1d histogram and 2d histogram to be user-defined intervals
			Combine the interval with the variable names in onw file.
            Separate 1d and 2d histogram interval in 2 files with all variables.

V6 Updates: Add the flexible input of sampling rate, polygon region and grid sizes of Lat & Lon

V7 Updates: Change the sampling rate starts from 3 and 4 count from 1 (Here we count from 0 so it starts from 2 and 3).
			Add the flexible input of data path, date and time period for averaging

V8 Updates: Change histogram count from averaged value to be all pixel values
			Added scale_factor, add_offset and _Fillvalue as attributes to each variabes.
			Found the problem in reading temperature and fixed it by changing the reading way:
				For netCDF4, the variable is done by (rdval * scale) + offst
				For MODIS HDF4 file, the variable should be done by (rdval-offst)*scale
				It needs to be reverted from the netCDF4 reading first, then convert it in the way of HDF file.
"""

import os
import sys
import h5py
import timeit
import random
import numpy as np
import pandas as pd
#from mpi4py import MPI
from netCDF4 import Dataset
from collections import OrderedDict
from datetime import date, datetime
from dateutil.rrule import rrule, DAILY, MONTHLY
from MODIS_Aggregation import *

if __name__ =='__main__':
# This is the main program for using concurrent to speed up the whole process

	#--------------STEP 1: Read User Inputs and Initial Paramters for Aggregation--------------------
	fname1,fname2,day_in_year,shift_hour,NTA_lats,NTA_lons,map_lon,map_lat,grid_lon,grid_lat,gap_x,gap_y,filenum, \
	grid_data,sts_switch,varnames,intervals_1d,intervals_2d,bin_num1,bin_num2,var_idx,spl_num,sts_name,histnames, \
	output_dir,l3name,unit_list,scale_list,offst_list,longname_list,fillvalue_list = read_user_inputs()

	#--------------STEP 2: Start Internal Aggregation------------------------------------------------

	# Start counting operation time
	start_time = timeit.default_timer() 

	#print(fname1.shape,fname2.shape)
	#for i in range(len(fname1)):
	#	print(fname1[i],fname2[i])
	#sys.exit()

	grid_data = run_modis_aggre(fname1,fname2,day_in_year,shift_hour,NTA_lats,NTA_lons,grid_lon,grid_lat,gap_x,gap_y,filenum, \
								grid_data,sts_switch,varnames,intervals_1d,intervals_2d,var_idx, spl_num, sts_name, histnames)
		
	# Compute the mean cloud fraction & Statistics (Include Min & Max & Standard deviation)

	# Reference for statstic parameters
	# sts_name[0]: min
	# sts_name[1]: max
	# sts_name[2]: mean -> total_value / count
	# sts_name[3]: count
	# sts_name[4]: square -> ((total_square_value / count) - mean**2) ** 0.5
	# sts_name[5]: histogram
	# sts_name[6]: joint histogram

	sts_idx = np.array(np.where(sts_switch == True))[0]
	print("Index of User-defined Statistics:",sts_idx)
	key_idx = 0
	for key in varnames:
		for i in sts_idx:
			if i == 0:
				grid_data[key+'_'+sts_name[0]] = grid_data[key+'_'+sts_name[0]].reshape([grid_lat,grid_lon])
			elif i == 1:
				grid_data[key+'_'+sts_name[1]] = grid_data[key+'_'+sts_name[1]].reshape([grid_lat,grid_lon])
			elif i == 2:
				grid_data[key+'_'+sts_name[2]] = (grid_data[key+'_'+sts_name[2]] / grid_data[key+'_'+sts_name[3]])
				grid_data[key+'_'+sts_name[2]] =  grid_data[key+'_'+sts_name[2]].reshape([grid_lat,grid_lon])
			elif i == 3:
				grid_data[key+'_'+sts_name[3]] =  grid_data[key+'_'+sts_name[3]].reshape([grid_lat,grid_lon])
			elif i == 4:
				grid_data[key+'_'+sts_name[4]] = ((grid_data[key+'_'+sts_name[4]] / grid_data[key+'_'+sts_name[3]].ravel()) - grid_data[key+'_'+sts_name[2]].ravel()**2)**0.5
				grid_data[key+'_'+sts_name[4]] = grid_data[key+'_'+sts_name[4]].reshape([grid_lat,grid_lon])
			elif i == 5:
				grid_data[key+'_'+sts_name[5]] = grid_data[key+'_'+sts_name[5]].reshape([grid_lat,grid_lon,bin_num1[key_idx]])
			elif i == 6:
				grid_data[key+'_'+sts_name[6]+histnames[key_idx]] = grid_data[key+'_'+sts_name[6]+histnames[key_idx]].reshape([grid_lat,grid_lon,bin_num1[key_idx],bin_num2[key_idx]])
		key_idx += 1	

	end_time = timeit.default_timer()

	
	print ("Operation Time in {:7.2f} seconds".format(end_time - start_time))
	
	#--------------STEP 3:  Create HDF5 file to store the result------------------------------
	ff=h5py.File(output_dir+l3name,'w')

	PC=ff.create_dataset('lat_bnd',data=map_lat)
	PC.attrs['units']='degrees'
	PC.attrs['long_name']='Latitude_boundaries'    

	PC=ff.create_dataset('lon_bnd',data=map_lon)
	PC.attrs['units']='degrees'
	PC.attrs['long_name']='Longitude_boundaries'    

	PCentry=ff.create_dataset('GRID_Counts',data=grid_data['GRID_Counts'].reshape([grid_lat,grid_lon]))
	PCentry.dims[0].label='lat_bnd'
	PCentry.dims[1].label='lon_bnd'
	PC.attrs['units']='none'
	PC.attrs['long_name']='grid_point_counts'    

	for i in range(sts_idx.shape[0]):
		cnt = 0
		for key in grid_data:

			if key.find("1km") != -1: 
				new_name = key.replace("_1km", "")
			else: 
				new_name = key

			if (sts_name[sts_idx[i]] in key) == True:  
				#print(sts_name[sts_idx[i]],key,grid_data[key].shape)
				#print(longname_list[cnt][:20],new_name)
				if sts_idx[i] >= 5:
					#print(cnt,intervals_1d[cnt],intervals_2d[cnt])
					addGridEntry(ff,new_name,unit_list[cnt],longname_list[cnt],fillvalue_list[cnt],scale_list[cnt],offst_list[cnt],grid_data[key],intervals_1d[cnt],intervals_2d[cnt])
				else:
					addGridEntry(ff,new_name,unit_list[cnt],longname_list[cnt],fillvalue_list[cnt],scale_list[cnt],offst_list[cnt],grid_data[key],intervals_1d[0],intervals_2d[0])
				cnt += 1
	
	ff.close()

	print(l3name+' Saved!')
	#---------------------------COMPLETED------------------------------------------------------

