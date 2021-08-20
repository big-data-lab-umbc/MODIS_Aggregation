#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Main Program: Run MODIS AGGREGATION IN MPI WITH FLEXIBLE STATISTICS

Created on 2020

@author: Jianyu Zheng (Email: jzheng3@umbc.edu)
"""

import os
import sys
import h5py
import timeit
import random
import calendar
import numpy as np
import pandas as pd
from mpi4py import MPI
from netCDF4 import Dataset
from collections import OrderedDict
from datetime import date, datetime
from dateutil.rrule import rrule, DAILY, MONTHLY
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

	print("-------- START AGGREGATION --------")
	
	# Initiate MPI
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	random.seed(rank)

	# Initiate the number of files for MPI
	remain   = size-total_file_num%size

	files_part1 = np.arange(total_file_num + remain)
	tasks_part1 = np.array(np.split(files_part1,size))

	files_part2 = np.arange(total_file_num - tasks_part1[rank].size * (size-remain)) + tasks_part1[rank].size * (size-remain)
	tasks_part2 = np.array(np.split(files_part2,remain))

	if rank < (size-remain):
		fileloop = tasks_part1[rank]
	else:
		fileloop = tasks_part2[rank-(size-remain)]

	print("Process {} calculating files from {} to {}... (Total: {} / {})".format(rank, fileloop[0],fileloop[-1],fileloop.shape[0],total_file_num))

	if rank == 0:
		grid_data = run_modis_aggre(fname1,fname2,day_in_year,shift_hour,NTA_lats,NTA_lons,grid_lon,grid_lat,gap_x,gap_y,fileloop, \
								    grid_data,sts_switch,varnames,intervals_1d,intervals_2d,var_idx,spl_num,sts_name,histnames)

		for i in range(1,size):
			results = comm.recv(source=i, tag=0)
			grid_data = addCounter(grid_data, results)

		# Compute the mean cloud fraction & Statistics (Include Min & Max & Standard deviation)

		# Reference for statstic parameters
		# sts_name[0]: min
		# sts_name[1]: max
		# sts_name[2]: mean / total_value
		# sts_name[3]: count
		# sts_name[4]: square
		# sts_name[5]: histogram
		# sts_name[6]: joint histogram

		sts_idx = np.array(np.where(sts_switch == True))[0]
		print("Index of User-defined Statistics:",sts_idx)
		print(grid_data['GRID_Counts'].reshape([grid_lat,grid_lon]))
		key_idx = 0
		for key in varnames:
			for i in sts_idx:
				if i == 0:
					grid_data[key+'_'+sts_name[0]] = grid_data[key+'_'+sts_name[0]].reshape([grid_lat,grid_lon])
				elif i == 1:
					grid_data[key+'_'+sts_name[1]] = grid_data[key+'_'+sts_name[1]].reshape([grid_lat,grid_lon])
				elif i == 2:
					grid_data[key+'_'+sts_name[2]] = grid_data[key+'_'+sts_name[2]] / grid_data[key+'_'+sts_name[3]]
					grid_data[key+'_'+sts_name[2]] = grid_data[key+'_'+sts_name[2]].reshape([grid_lat,grid_lon])
				elif i == 3:
					grid_data[key+'_'+sts_name[3]] = grid_data[key+'_'+sts_name[3]].reshape([grid_lat,grid_lon])
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
		ff=h5py.File(output_dir+l3name+'MPI','w')

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
					if sts_idx[i] >= 5:
						addGridEntry(ff,new_name,unit_list[cnt],longname_list[cnt],fillvalue_list[cnt],scale_list[cnt],offst_list[cnt],grid_data[key],intervals_1d[cnt],intervals_2d[cnt])
					else:
						addGridEntry(ff,new_name,unit_list[cnt],longname_list[cnt],fillvalue_list[cnt],scale_list[cnt],offst_list[cnt],grid_data[key],intervals_1d[0],intervals_2d[0])
					cnt += 1

		ff.close()

		print(l3name+' Saved!')
		print("-------- AGGREGATION COMPLETED --------")

	else:
		results = run_modis_aggre(fname1,fname2,day_in_year,shift_hour,NTA_lats,NTA_lons,grid_lon,grid_lat,gap_x,gap_y,fileloop, \
								  grid_data,sts_switch,sts_name,histnames,varnames,intervals_1d,intervals_2d,var_idx,spl_num)
		massage = "Process {} finished".format(rank)
		print(massage)
		comm.send(results, dest=0, tag=0)


#---------------------------COMPLETED------------------------------------------------------
