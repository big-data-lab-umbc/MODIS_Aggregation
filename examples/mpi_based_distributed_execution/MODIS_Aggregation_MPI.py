#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Main Program: Run MODIS AGGREGATION IN MPI WITH FLEXIBLE STATISTICS
Created on 2019
@author: Jianyu Zheng
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

	#-------------STEP 0: Read the input from User --------
	# checking user input
	if (len(sys.argv) != 16) & (len(sys.argv) != 17):
		print("## ERROR!!!")
		print("## Wrong user input")
		print("## Please see the example below:")
		print("## python MODIS_Aggregation_Serial.py \n" + \
			  "   <Data Path> <Start Date: yyyy-mm-dd> <End Date: yyyy-mm-dd> \n" + \
			  "   <Polygon boundaries> <Lat & Lon Grid Size > \n" + \
			  "   <Sampling number larger than 0> \n" + \
			  "   <1/0> <1/0> <1/0> \n" + \
			  "   <1/0> <1/0> <1/0> \n" + \
			  "   <1/0> <Variable Imput File> <JHist Variable Imput File>")

		#start_date=np.fromstring(sys.argv[2], dtype=np.int, sep='/' )
		#end_date=np.fromstring(sys.argv[2], dtype=np.int, sep='/' )
		#print("Date:",start_date[0],start_date[1],start_date[2])
		#spl_num = np.int(sys.argv[3][1:-1])
		#poly=np.fromstring(sys.argv[1][1:-1], dtype=np.int, sep=',' )
		#grid=np.fromstring(sys.argv[2][1:-1], dtype=np.float, sep=',' )
		#print(spl_num,poly,grid)
		sys.exit()
	else:
		# Define the sampling rate, boundaries of the selected polygon region & the grid size of Lat & Lon
		spl_num = np.int(sys.argv[6][1:-1])
		poly=np.fromstring(sys.argv[4][1:-1], dtype=np.int, sep=',' )
		grid=np.fromstring(sys.argv[5][1:-1], dtype=np.float, sep=',' )

		# Define the statistics names for HDF5 output
		sts_name = ['Minimum','Maximum','Mean','Pixel_Counts', \
					'Standard_Deviation','Histogram_Counts','Jhisto_vs_']

		# Pass system arguments to the function
		sts_switch = np.array(sys.argv[7:14],dtype=np.int)
		sts_switch = np.array((sts_switch == 1))
		varlist  = sys.argv[14]

		# Read the variable names from the variable name list
		text_file = np.array(pd.read_csv(varlist, header=0, delim_whitespace=True)) #open(varlist, "r")
		varnames  = text_file[:,0]

		if sts_switch[5] == True:
			intervals_1d = text_file[:,1] # This is a string interval arrays
		else:
			intervals_1d = [0]

		if sts_switch[6] == True:
			# Read the joint histogram names from the variable name list
			jvarlist = sys.argv[15]
			text_file = np.array(pd.read_csv(jvarlist, header=0, delim_whitespace=True)) #open(varlist, "r")
			histnames = text_file[:,1]
			var_idx   = text_file[:,2] #This is the index of the input variable name which is used for 2D histogram
			intervals_2d = text_file[:,3]
		else:
			intervals_2d,var_idx,histnames = [0],[0],np.nan

	#-------------STEP 1: Set up the specific directory --------
	data_path_file = np.array(pd.read_csv(sys.argv[1], header=0, delim_whitespace=True))
	MYD06_dir    = data_path_file[0,0] #'/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD06_L2/'
	MYD06_prefix = data_path_file[0,1] #'MYD06_L2.A'
	MYD03_dir    = data_path_file[1,0] #'/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD03/'
	MYD03_prefix = data_path_file[1,1] #'MYD03.A'
	fileformat = 'hdf'

	output_path_file = np.array(pd.read_csv(sys.argv[1], header=3, delim_whitespace=True))
	output_dir = output_path_file[0,0]
	output_prefix = output_path_file[0,1]

	#-------------STEP 2: Set up spactial and temporal resolution & variable names----------
	NTA_lats = [poly[0],poly[1]] #[  0,40] #[-90,90]   #[-30,30]
	NTA_lons = [poly[2],poly[3]] #[-40,60] #[-180,180] #[-60,60]

	gap_x, gap_y = grid[1],grid[0] #0.5,0.625

	if ((NTA_lons[-1]-NTA_lons[0])%gap_x != 0) | ((NTA_lats[-1]-NTA_lats[0])%gap_y != 0):
		print("## ERROR!!!")
		print("## Grid size should be dividable by the dimension of the selected region.")
		print("## If you choose the region of latitude  from -40 to 40, then you gird size (Latitude ) should be dividable by 80.")
		print("## If you choose the region of longitude from  20 to 35, then you gird size (Longitude) should be dividable by 55.")
		print("## Please try again!")
		sys.exit()

	map_lon = np.arange(NTA_lons[0],NTA_lons[1],gap_x)
	map_lat = np.arange(NTA_lats[0],NTA_lats[1],gap_y)
	Lon,Lat = np.meshgrid(map_lon,map_lat)
	grid_lon=np.int((NTA_lons[-1]-NTA_lons[0])/gap_x)
	grid_lat=np.int((NTA_lats[-1]-NTA_lats[0])/gap_y)

	#--------------STEP 3: Create arrays for level-3 statistics data-------------------------
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

	#--------------STEP 4: Read the filename list for different time period-------------------
	fname1,fname2 = [],[]

	start_date = np.fromstring(sys.argv[2], dtype=np.int, sep='/' )
	end_date   = np.fromstring(sys.argv[3], dtype=np.int, sep='/' )
	start = date(start_date[0], start_date[1], start_date[2])
	until = date(end_date[0], end_date[1], end_date[2])

	for dt in rrule(DAILY, interval=1, dtstart=start, until=until):
		year  = np.array([np.int(dt.strftime("%Y"))])
		month = np.array([np.int(dt.strftime("%m"))])
		day   = np.array([np.int(dt.strftime("%d"))])
		time  = np.arange(24) #np.int(dt.strftime("%H"))

		daynew = dt.toordinal()
		yearstart = datetime(year,1,1)
		yearend   = calendar.monthrange(year, 12)[1]

		day_yearstart = yearstart.toordinal()
		day_yearend = datetime(year,12,yearend).toordinal()

		day_in_year = np.array([(daynew-day_yearstart)+1])
		end_in_year = np.array([(day_yearend-day_yearstart)+1])

		# Adjust to 3 hours previous/after the End Date for the orbit gap/overlap problem
		if (dt.year == until.year) & (dt.month == until.month) & (dt.day == until.day):
			shift_hour = 3
			time    = np.append(np.arange(24),np.arange(shift_hour))
			year    = [year[0],year[0]]
			day_in_year = [day_in_year[0],day_in_year[0] + 1]
			if day_in_year[1] > end_in_year:
				year[1]   -= 1
				yearstart = datetime(year[1],1,1)
				yearend   = datetime(year[1],12,31)
				day_yearstart = yearstart.toordinal()
				day_yearend   = yearend.toordinal()
				day_in_year[1] = (day_yearend-day_yearstart)+1

		# Start reading Level-2 files
		fname_tmp1,fname_tmp2 = read_filelist(MYD06_dir,MYD06_prefix,MYD03_dir,MYD03_prefix,year,day_in_year,time,fileformat)
		fname1 = np.append(fname1,fname_tmp1)
		fname2 = np.append(fname2,fname_tmp2)

		#print(fname1.shape,fname2.shape)

	total_file_num = len(fname1) #np.arange(len(fname1))
	#print(fname1)

	#--------------STEP 5: Read Attributes of each variables----------------------------------
	unit_list = []
	scale_list = []
	offst_list = []
	longname_list = []
	fillvalue_list = []

	ncfile=Dataset(fname1[0],'r')

	# Read the User-defined variables from MYD06 product
	tmp_idx = 0
	cloud_fraction_flag = 0
	for key in varnames:
		if key == 'cloud_fraction':
			name_idx = tmp_idx
			cloud_fraction_flag = 1
			continue #Ignoreing Cloud_Fraction from the input file
		else:
			tmp_data,data_dim,lonam,unit,fill,scale,offst = readEntry(key,ncfile,spl_num)
			unit_list  = np.append(unit_list,unit)
			scale_list = np.append(scale_list,scale)
			offst_list = np.append(offst_list,offst)
			longname_list = np.append(longname_list, lonam)
			fillvalue_list = np.append(fillvalue_list, fill)
			tmp_idx += 1

	if cloud_fraction_flag == 0:
		print("## ERROR!!!")
		print("## The cloud fraction varibale is not included in the aggregation.")
		print("## Please add the cloud fraction to the input_file.csv ")
		sys.exit()

	# Add the long name of cloud freaction at the first row
	CM_unit     = 'none'
	CM_longname = 'Cloud Fraction from Cloud Mask (cloudy & prob cloudy)'
	CM_fillvalue = -9999
	CM_scale_factor = 0.0001
	CM_add_offset   = 0.0
	unit_list      = np.insert(unit_list,      name_idx, CM_unit)
	scale_list     = np.insert(scale_list,     name_idx, CM_scale_factor)
	offst_list     = np.insert(offst_list,     name_idx, CM_add_offset)
	longname_list  = np.insert(longname_list,  name_idx, CM_longname)
	fillvalue_list = np.insert(fillvalue_list, name_idx, CM_fillvalue)

	ncfile.close()

	#--------------STEP 6: Start Aggregation------------------------------------------------

	# Start counting operation time
	start_time = timeit.default_timer()

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
		#print(rank)
		fileloop = tasks_part1[rank]
	else:
		#print(rank,rank-(size-remain))
		fileloop = tasks_part2[rank-(size-remain)]

	print("Process {} calculating files from {} to {}... (Total: {} / {})".format(rank, fileloop[0],fileloop[-1],fileloop.shape[0],total_file_num))

	#sys.exit()

	if rank == 0:
		#grid_data = run_modis_aggre(fname1,fname2,day_in_year,shift_hour,NTA_lats,NTA_lons,grid_lon,grid_lat,gap_x,gap_y,fileloop, \
		#							grid_data,sts_switch,varnames,intervals_1d,intervals_2d,var_idx, spl_num, sts_name, histnames)
		grid_data = run_modis_aggre(fname1,fname2,day_in_year,shift_hour,NTA_lats,NTA_lons,grid_lon,grid_lat,gap_x,gap_y,fileloop, \
								    grid_data,sts_switch,varnames,intervals_1d,intervals_2d,var_idx, spl_num, sts_name  ,histnames)

		for i in range(1,size):
			results = comm.recv(source=i, tag=0)
			grid_data = addCounter(grid_data, results)

		#counterSumOp = MPI.Op.Create(addCounter, commute=True)
		#grid_data = comm.allreduce(results, op=counterSumOp)

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

		#print('Mean_Fraction:')
		#print( Mean_Fraction  )

		print ("Operation Time in {:7.2f} seconds".format(end_time - start_time))

		#--------------STEP 7:  Create HDF5 file to store the result------------------------------
		l3name  = output_prefix + '.A{:04d}{:03d}.'.format(year[0],day_in_year[0])

		subname = 'MPI_output_monthly_5km.h5'
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
					#print(cnt,sts_name[sts_idx[i]],key,grid_data[key].shape)
					#print(longname_list[cnt][:20],new_name)
					if sts_idx[i] >= 5:
						#print(cnt,unit_list.shape,intervals_1d[cnt],intervals_2d[cnt])
						addGridEntry(ff,new_name,unit_list[cnt],longname_list[cnt],fillvalue_list[cnt],scale_list[cnt],offst_list[cnt],grid_data[key],intervals_1d[cnt],intervals_2d[cnt])
					else:
						#print(cnt,unit_list.shape,intervals_1d[cnt],intervals_2d[cnt])
						addGridEntry(ff,new_name,unit_list[cnt],longname_list[cnt],fillvalue_list[cnt],scale_list[cnt],offst_list[cnt],grid_data[key],intervals_1d[0],intervals_2d[0])
					cnt += 1

		ff.close()

		print(l3name+subname+' Saved!')

	else:
		results = run_modis_aggre(fname1,fname2,day_in_year,shift_hour,NTA_lats,NTA_lons,grid_lon,grid_lat,gap_x,gap_y,fileloop, \
								  grid_data,sts_switch,varnames,intervals_1d,intervals_2d,var_idx, spl_num, sts_name  ,histnames)
		massage = "Process {} finished".format(rank)
		print(massage)
		comm.send(results, dest=0, tag=0)


#---------------------------COMPLETED------------------------------------------------------
