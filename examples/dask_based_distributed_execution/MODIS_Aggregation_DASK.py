#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Main Program: Run MODIS AGGREGATION IN DASK WITH FLEXIBLE STATISTICS
Created on 2021
@author: Xin Huang
"""

import os
import sys
import timeit
import random
import calendar
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import h5py
from collections import OrderedDict
from datetime import date, datetime
from dateutil.rrule import rrule, DAILY, MONTHLY
import dask
from dask.distributed import as_completed
from dask_jobqueue import SLURMCluster
from dask.distributed import Client, LocalCluster
from dask.distributed import wait
from MODIS_Aggregation import *


if __name__ =='__main__':
# This is the main program for using concurrent to speed up the whole process

	#--------------STEP 1: Read User Inputs and Initial Paramters for Aggregation--------------------
	fname1,fname2,day_in_year,shift_hour,NTA_lats,NTA_lons,map_lon,map_lat,grid_lon,grid_lat,gap_x,gap_y,total_file_num, \
	grid_data,sts_switch,varnames,intervals_1d,intervals_2d,bin_num1,bin_num2,var_idx,spl_num,sts_name,histnames, \
	output_dir,l3name,unit_list,scale_list,offst_list,longname_list,fillvalue_list = read_user_inputs()

	#--------------STEP 2: Start Aggregation------------------------------------------------

	# Start counting operation time
	start_time = timeit.default_timer()

	#print(fname1.shape,fname2.shape)
	#for i in range(len(fname1)):
	#	print(fname1[i],fname2[i])
	#sys.exit()

	# parallel implement
	file_num = len(fname1)
	CHUNK_NUM = 256
	print("file_num:", file_num)
	#files_range = range(file_num)
	#chunks = np.array_split(files_range, CHUNK_NUM)
	#print("chunks:", chunks)
	print("filename1:", fname1)
	print("filename2:", fname2)
	filename1_chunks = np.array_split(fname1, CHUNK_NUM)
	filename2_chunks = np.array_split(fname2, CHUNK_NUM)

	# Initiate and process the parallel by Dask

	kwargv = { "day_in_year": day_in_year, "shift_hour": shift_hour, "NTA_lats": NTA_lats, "NTA_lons": NTA_lons, "grid_lon": grid_lon,"grid_lat": grid_lat, "gap_x": gap_x, "gap_y": gap_y, "filenum": filenum, "sts_switch":sts_switch, "sts_name":sts_name, "histnames":histnames, "varnames": varnames, "intervals_1d":intervals_1d, "intervals_2d":intervals_2d, "var_idx":var_idx, "spl_num": spl_num}

	cluster = SLURMCluster(cores=32, memory='390 GB',processes=32, project='pi_jianwu',\
		queue='high_mem', walltime='16:00:00', job_extra=['--exclusive', '--qos=medium+'])
	cluster.scale(jobs=8)
	print(cluster.job_script())
	client = Client(cluster)
	#print("varnames:", varnames)
	#print("")
	#print("kwargv:", kwargv)
	#client = Client()
	tt = client.map(run_modis_aggre, filename1_chunks, filename2_chunks, **kwargv)

	#init the global data result
	grid_data = {}
	bin_num1 = np.zeros(len(varnames)).astype(np.int)
	bin_num2 = np.zeros(len(varnames)).astype(np.int)
	grid_data['GRID_Counts'] = np.zeros(grid_lat*grid_lon).astype(np.int)
	key_idx = 0
	for key in varnames:
		if sts_switch[0] == True:
			grid_data[key+'_'+sts_name[0]] = np.zeros(grid_lat*grid_lon) + np.inf
		if sts_switch[1] == True:
			grid_data[key+'_'+sts_name[1]] = np.zeros(grid_lat*grid_lon) - np.inf
		if (sts_switch[2] == True) | (sts_switch[3] == True) | (sts_switch[4] == True):
			grid_data[key+'_'+sts_name[2]] = np.zeros(grid_lat*grid_lon)
			grid_data[key+'_'+sts_name[3]] = np.zeros(grid_lat*grid_lon)
			grid_data[key+'_'+sts_name[4]] = np.zeros(grid_lat*grid_lon)
		if sts_switch[5] == True:
			bin_interval1 = np.fromstring(intervals_1d[key_idx], dtype=np.float, sep=',' )
			bin_num1[key_idx] = bin_interval1.shape[0]-1
			grid_data[key+'_'+sts_name[5]] = np.zeros((grid_lat*grid_lon,bin_num1[key_idx]))

			if sts_switch[6] == True:
				bin_interval2 = np.fromstring(intervals_2d[key_idx], dtype=np.float, sep=',' )
				bin_num2[key_idx] = bin_interval2.shape[0]-1
				grid_data[key+'_'+sts_name[6]+histnames[key_idx]] = np.zeros((grid_lat*grid_lon,bin_num1[key_idx],bin_num2[key_idx]))

		key_idx += 1

	# aggregate the result
	for future, result in as_completed(tt, with_results= True):
		for key in result:
			if key.find("Minimum") != -1:
				grid_data[key] = np.fmin(grid_data[key],result[key])
			elif key.find("Maximum") != -1:
				grid_data[key] = np.fmax(grid_data[key],result[key])
			else:
				grid_data[key] += result[key]


	#grid_data = run_modis_aggre(fname1,fname2,day_in_year,shift_hour,NTA_lats,NTA_lons,grid_lon,grid_lat,gap_x,gap_y,filenum, \
	#							grid_data,sts_switch,varnames,intervals_1d,intervals_2d,var_idx, spl_num, sts_name, histnames)

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
				if key == 'cloud_fraction':
					grid_data[key+'_'+sts_name[4]] = ((grid_data[key+'_'+sts_name[4]] / grid_data['GRID_Counts']) - grid_data[key+'_'+sts_name[2]].ravel()**2)**0.5
				else:
					grid_data[key+'_'+sts_name[4]] = grid_data[key+'_'+sts_name[4]] / grid_data['GRID_Counts']
				grid_data[key+'_'+sts_name[4]] = grid_data[key+'_'+sts_name[4]].reshape([grid_lat,grid_lon])
			elif i == 5:
				grid_data[key+'_'+sts_name[5]] = grid_data[key+'_'+sts_name[5]].reshape([grid_lat,grid_lon,bin_num1[key_idx]])
			elif i == 6:
				grid_data[key+'_'+sts_name[6]+histnames[key_idx]] = grid_data[key+'_'+sts_name[6]+histnames[key_idx]].reshape([grid_lat,grid_lon,bin_num1[key_idx],bin_num2[key_idx]])
		key_idx += 1

	end_time = timeit.default_timer()

	#print('Mean_Fraction:')
	#print( Mean_Fraction  )

	print ("Operation Time in {:7.2f} seconds".format(end_time - start_time))

	#--------------STEP 7:  Create HDF5 file to store the result------------------------------
	subname = 'dask_output_1month_chunk_n8_chunk256.h5'
	ff=h5py.File(output_dir+l3name+subname,'w')

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

	print(l3name+subname+' Saved!')
	#---------------------------COMPLETED------------------------------------------------------
