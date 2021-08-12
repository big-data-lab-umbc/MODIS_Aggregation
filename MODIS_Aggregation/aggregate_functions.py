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
from mpi4py import MPI
from netCDF4 import Dataset
from collections import OrderedDict
from datetime import date, datetime
from dateutil.rrule import rrule, DAILY, MONTHLY
	
def read_filelist(loc_dir,prefix,yr,day,fileformat):
	# Read the filelist in the specific directory
	str = os.popen("ls "+ loc_dir + prefix + yr + day + "*."+fileformat).read()
	fname = np.array(str.split("\n"))
	fname = np.delete(fname,len(fname)-1)

	return fname

def readEntry(key,ncf):
	# Read the MODIS variables based on User's name list
	rdval=np.array(ncf.variables[key]).astype(np.float)

	#For netCDF4, the variable is done by (rdval * scale) + offst
	#For MODIS HDF4 file, the variable should be done by (rdval-offst)*scale 
	#It needs to be reverted from the netCDF4 reading first, then convert it in the way of HDF file.

	# Read the attributes of the variable
	unit = ncf.variables[key].units
	scale= ncf.variables[key].scale_factor
	offst= ncf.variables[key].add_offset
	lonam= ncf.variables[key].long_name
	fillvalue = ncf.variables[key]._FillValue

	rdval[np.where(rdval == fillvalue)] = np.nan

	#Sampling the variable
	rdval = rdval[2::spl_num,3::spl_num]

	return rdval,lonam,unit,fillvalue,scale,offst

def read_MODIS(varnames,fname1,fname2): 
	# Store the data from variables after reading MODIS files
	data={}
	
	# Read the Cloud Mask from MYD06 product
	ncfile=Dataset(fname1,'r')

	#CM1km = readEntry('Cloud_Mask_1km',ncfile)
	#CM1km = np.array(ncfile.variables['Cloud_Mask_1km'])
	#data['CM'] = (np.array(CM1km[:,:,0],dtype='byte') & 0b00000110) >>1
	
	d06_CM = ncfile.variables['Cloud_Mask_1km'][:,:,0]
	CM1km  = d06_CM[2::spl_num,3::spl_num]
	data['CM'] = (np.array(CM1km,dtype='byte') & 0b00000110) >>1
	data['CM'] = data['CM'].astype(np.float)

	# Read the User-defined variables from MYD06 product
	for key in varnames:
		if key == 'cloud_fraction': 
			continue #Ignoreing Cloud_Fraction from the input file
		else:
			data[key],lonam,unit,fill,scale,offst = readEntry(key,ncfile)
			data[key] = (data[key] - offst) / scale 
			data[key] = (data[key] - offst) * scale 
			
	ncfile.close()

	# Read the common variables (Latitude & Longitude) from MYD03 product
	ncfile=Dataset(fname2,'r')
	d03_lat  = np.array(ncfile.variables['Latitude'][:,:])
	d03_lon  = np.array(ncfile.variables['Longitude'][:,:])
	lat  = d03_lat[2::spl_num,3::spl_num]
	lon  = d03_lon[2::spl_num,3::spl_num]
	attr_lat = ncfile.variables['Latitude']._FillValue
	attr_lon = ncfile.variables['Longitude']._FillValue

	# If the variable is not 1km product, exit and tell the User to reset the variables.
	for key in varnames:
		if key == 'cloud_fraction': continue #Ignoreing Cloud_Fraction from the input file	
		if data[key].shape[0] != lat.shape[0]:
			print("The dimension of varibale '"+key+"' is not match with latitude & longitude.")
			print("Input variables should have 1km resolution.")
			print("Check your varibales.")
			sys.exit()

	#Use _FillValue to remove fill data in lat & lon
	lat[np.where(lat == attr_lat)] = np.nan
	lon[np.where(lat == attr_lat)] = np.nan
	data['CM'][np.where(lat == attr_lat)] = np.nan #which will not be identified by lines 80-83 

	lat[np.where(lon == attr_lon)] = np.nan
	lon[np.where(lon == attr_lon)] = np.nan
	data['CM'][np.where(lon == attr_lon)] = np.nan #which will not be identified by lines 80-83
	ncfile.close()

	return lat,lon,data

def cal_stats(z,key,grid_data,min_val,max_val,tot_val,count,all_val,all_val_2d, \
			  sts_switch,sts_name,intervals_1d,intervals_2d,key_idx):
# Calculate Statistics pamameters
					
	#Min and Max
	if sts_switch[0] == True:
		if  grid_data[key+'_'+sts_name[0]][z] > min_val:
			grid_data[key+'_'+sts_name[0]][z] = min_val
	
	if sts_switch[1] == True:
		if  grid_data[key+'_'+sts_name[1]][z] < max_val:
			grid_data[key+'_'+sts_name[1]][z] = max_val

	#Total and Count for Mean
	if (sts_switch[2] == True) | (sts_switch[3] == True):
		grid_data[key+'_'+sts_name[2]][z] += tot_val
		grid_data[key+'_'+sts_name[3]][z] += count
	
	#Standard Deviation 
	if sts_switch[4] == True:
		grid_data[key+'_'+sts_name[4]][z] += tot_val**2

	#1D Histogram 
	if sts_switch[5] == True:	
		bin_interval1 = np.fromstring(intervals_1d[key_idx], dtype=np.float, sep=',' )
		if all_val.size == 1: 
			all_val = np.array([all_val])
		else:
			hist_idx1 = np.histogram(all_val,bins=bin_interval1)[0]
			grid_data[key+'_'+sts_name[5]][z,:] += hist_idx1
		#for i in range(all_val.size):
		#	hist_idx1 = np.where(bin_interval1 <= all_val[i])[0]
		#	#hist_idx1 = 0 if len(hist_idx1) == 0 else hist_idx1[-1]
		#	#if hist_idx1 > (grid_data[key+'_'+sts_name[5]].shape[1]-1): hist_idx1 = (grid_data[key+'_'+sts_name[5]].shape[1]-1) 
		#	grid_data[key+'_'+sts_name[5]][z, hist_idx1] += 1
		
	#2D Histogram 
	if sts_switch[6] == True:
		bin_interval1 = np.fromstring(intervals_1d[key_idx], dtype=np.float, sep=',' )
		bin_interval2 = np.fromstring(intervals_2d[key_idx], dtype=np.float, sep=',' )
		if all_val.size == 1:
			all_val = np.array([all_val])
			all_val_2d = np.array([all_val_2d])
		else:
			hist_idx2 = np.histogram2d(all_val,all_val_2d,bins=(bin_interval1,bin_interval2))[0]
			grid_data[key+'_'+sts_name[6]+histnames[key_idx]][z,:,:] += hist_idx2

		#for i in range(all_val_2d.size):
		#	hist_idx1 = np.where(bin_interval1 <= all_val[i])[0]
		#	hist_idx1 = 0 if len(hist_idx1) == 0 else hist_idx1[-1]
		#	if hist_idx1 > (grid_data[key+'_'+sts_name[5]].shape[1]-1): hist_idx1 = (grid_data[key+'_'+sts_name[5]].shape[1]-1) 
		#
		#	hist_idx2 = np.where(bin_interval2 <= all_val_2d[i])[0]
		#	hist_idx2 = 0 if len(hist_idx2) == 0 else hist_idx2[-1]
		#	if hist_idx2 > (grid_data[key+'_'+sts_name[6]+histnames[key_idx]].shape[2]-1): hist_idx2 = (grid_data[key+'_'+sts_name[6]+histnames[key_idx]].shape[2]-1) 
		#
		#	grid_data[key+'_'+sts_name[6]+histnames[key_idx]][z, hist_idx1,hist_idx2] += 1
		
	return grid_data

def run_modis_aggre(fname1,fname2,NTA_lats,NTA_lons,grid_lon,grid_lat,gap_x,gap_y,hdfs, \
					grid_data,sts_switch,varnames,intervals_1d,intervals_2d,var_idx):
	# This function is the data aggregation loops by number of files
	hdfs = np.array(hdfs)
	for j in hdfs:#range(1):#hdfs:
		print("File Number: {} / {}".format(j,hdfs[-1]))
	
		# Read Level-2 MODIS data
		lat,lon,data = read_MODIS(varnames,fname1[j],fname2[j])
		CM = data['CM']
		
		# Restrain lat & lon & variables in the required region 
		res_idx = np.where((lat > NTA_lats[0]) & (lat < NTA_lats[1]) & (lon > NTA_lons[0]) & (lon < NTA_lons[1]))
		lat = lat[res_idx]
		lon = lon[res_idx]
		CM  = CM [res_idx]
		
		# Ravel the 2-D data to 1-D array
		lat = lat.ravel()
		lon = lon.ravel()
		CM  = CM.ravel()

		key_idx = 0
		for key in varnames:
			if key == 'cloud_fraction': 
				CF_key_idx = key_idx
				key_idx += 1
				continue #Ignoreing Cloud_Fraction from the input file	
			data[key] = data[key][res_idx].ravel()
			key_idx += 1
			
		# Locate the lat lon index into 3-Level frid box
		idx_lon = np.round((lon-NTA_lons[0])/gap_x).astype(int)
		idx_lat = np.round((lat-NTA_lats[0])/gap_y).astype(int)

		latlon_index=(idx_lat*grid_lon)+idx_lon

		latlon_index_unique = np.unique(latlon_index)

		#print(lon[0],idx_lon[0],lat[0],idx_lat[0])
		#print(latlon_index_unique.max(),grid_lat*grid_lon)
		
		for i in np.arange(latlon_index_unique.size):
		#-----loop through all the grid boxes ocupied by this granule------#
			z=latlon_index_unique[i]
			if((z >= 0) & (z < (grid_lat*grid_lon))):

				# For cloud fraction
				TOT_pix = np.sum(CM[np.where(latlon_index == z)]>=0).astype(float)
				CLD_pix = np.sum(CM[np.where(latlon_index == z)]<=1).astype(float)

				#local_data = CM[np.where(latlon_index == z)]
				#if local_data[np.where(np.isnan(local_data) == 0)].size == 0: 
				#	print('All NaN is Ture.')

				Fraction = CLD_pix / TOT_pix

				if len(intervals_2d) != 1:	
					pixel_data_2d = data[varnames[var_idx[CF_key_idx]]]
					ave_val_2d = np.nansum(pixel_data_2d[np.where(latlon_index == z)]).astype(float) / TOT_pix
				else:
					ave_val_2d = 0

				# Calculate Statistics pamameters
				grid_data = cal_stats(z,"cloud_fraction",grid_data, \
									  Fraction,Fraction,CLD_pix,TOT_pix,Fraction,ave_val_2d, \
									  sts_switch,sts_name,intervals_1d,intervals_2d,CF_key_idx)

				# For other variables
				key_idx = 0
				for key in varnames:
					if key == 'cloud_fraction': #Ignoreing Cloud_Fraction from the input file	
						key_idx += 1
						continue 
					pixel_data = data[key]

					tot_val = np.nansum(pixel_data[np.where(latlon_index == z)]).astype(float)
					#ave_val = tot_val / TOT_pix
					all_val = np.array(pixel_data[np.where(latlon_index == z)]).astype(float)
					max_val = np.nanmax(pixel_data[np.where(latlon_index == z)]).astype(float)
					min_val = np.nanmin(pixel_data[np.where(latlon_index == z)]).astype(float)

					#local_data = pixel_data[np.where(latlon_index == z)]
					#print(local_data.size,z)
					#print(z,tot_val,max_val,min_val,CLD_pix,TOT_pix)
					#if local_data[np.where(np.isnan(local_data) == 0)].size == 0: 
					#	print('All NaN is Ture.',key,tot_val,max_val,min_val,CLD_pix,TOT_pix)
					
					if len(intervals_2d) != 1:
						pixel_data_2d = data[varnames[var_idx[key_idx]]]
						all_val_2d = np.array(pixel_data_2d[np.where(latlon_index == z)]).astype(float)
					else:
						all_val_2d = 0

					# Calculate Statistics pamameters
					grid_data = cal_stats(z,key,grid_data, \
										  min_val,max_val,tot_val,CLD_pix,all_val,all_val_2d, \
										  sts_switch,sts_name,intervals_1d,intervals_2d,key_idx)

					key_idx += 1

	return grid_data


def addGridEntry(f,name,units,long_name,fillvalue,scale_factor,add_offset,data):
	'''
	f:h5py.File()
	-------------------------------------
	Ex.
	self.addGridEntry(f,'CF','Fraction','Cloud_Fraction',total_cloud_fraction)

	For netCDF4, the variable is done by (rdval * scale) + offst
	For MODIS HDF4 file, the variable should be done by (rdval-offst)*scale 
	It needs to be reverted from the netCDF4 reading first, then convert it in the way of HDF file.
	'''
	if (('Histogram_Counts' in name) == True) | (('Jhisto_vs_' in name) == True) | (('Pixel_Counts' in name) == True):
		original_data = data.astype(np.int)
	elif (('Maximum' in name) == True) | (('Minimum' in name) == True):
		tmp_data = data/scale_factor + add_offset
		tmp_data[np.where(np.isinf(tmp_data) == 1)]=fillvalue
		original_data = tmp_data.astype(np.int)
	else: 
		tmp_data = data/scale_factor + add_offset
		tmp_data[np.where(np.isnan(tmp_data) == 1)]=fillvalue
		original_data = tmp_data.astype(np.int)

	PCentry=f.create_dataset(name,data=original_data)
	PCentry.dims[0].label='lat_bnd'
	PCentry.dims[1].label='lon_bnd'
	PCentry.attrs['units']       = np.str(units)
	PCentry.attrs["long_name"]   = np.str(long_name)
	PCentry.attrs['_FillValue']  = fillvalue
	PCentry.attrs['scale_factor']= scale_factor
	PCentry.attrs['add_offset']  = add_offset

if __name__ =='__main__':
# This is the main program for using concurrent to speed up the whole process
	
	#-------------STEP 0: Read the input from User --------
	# checking user input
	if (len(sys.argv) != 16) & (len(sys.argv) != 17):
		print("Wrong user input")
		print("usage: python aggre_stats_mpi.py <Data Path> <Start Date> <End Date> \
												<Polygon boundaries> <Lat & Lon Grid Size > \
												<Sampling number larger than 0> \
												<1/0> <1/0> <1/0> \
												<1/0> <1/0> <1/0> \
												<1/0> <Variable Imput File> <JHist Variable Imput File>")
		
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
			intervals_2d,var_idx = [0],[0]
		
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
		print("Grid size should be dividable by the dimension of the selected region.")
		print("If you choose the region of latitude  from -40 to 40, then you gird size (Latitude ) should be dividable by 80.")
		print("If you choose the region of longitude from  20 to 35, then you gird size (Longitude) should be dividable by 55.")
		print("Please try again!")
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
		year  = np.int(dt.strftime("%Y"))
		month = np.int(dt.strftime("%m"))
		day   = np.int(dt.strftime("%d"))

		data = datetime(year,month,day)
		daynew = data.toordinal()
		yearstart = datetime(year,1,1)
		day_yearstart = yearstart.toordinal()
		day_in_year = (daynew-day_yearstart)+1

		yc ='%04i' % year
		dc ='%03i' % day_in_year

		fname_tmp1 = read_filelist(MYD06_dir,MYD06_prefix,yc,dc,fileformat)
		fname_tmp2 = read_filelist(MYD03_dir,MYD03_prefix,yc,dc,fileformat)
		fname1 = np.append(fname1,fname_tmp1)
		fname2 = np.append(fname2,fname_tmp2)
	print(year,month)

	filenum = np.arange(len(fname1))
	print(len(fname1))

	#--------------STEP 5: Read Attributes of each variables----------------------------------
	unit_list = []
	scale_list = []
	offst_list = []
	longname_list = []
	fillvalue_list = []

	ncfile=Dataset(fname1[0],'r')

	# Read the User-defined variables from MYD06 product
	tmp_idx = 0
	for key in varnames:
		if key == 'cloud_fraction': 
			name_idx = tmp_idx
			continue #Ignoreing Cloud_Fraction from the input file
		else:
			tmp_data,lonam,unit,fill,scale,offst = readEntry(key,ncfile)
			unit_list  = np.append(unit_list,unit)
			scale_list = np.append(scale_list,scale)
			offst_list = np.append(offst_list,offst)
			longname_list = np.append(longname_list, lonam)
			fillvalue_list = np.append(fillvalue_list, fill)
			tmp_idx += 1

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

	grid_data = run_modis_aggre(fname1,fname2,NTA_lats,NTA_lons,grid_lon,grid_lat,gap_x,gap_y,filenum, \
								grid_data,sts_switch,varnames,intervals_1d,intervals_2d,var_idx)
		
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
				grid_data[key+'_'+sts_name[4]] =  grid_data[key+'_'+sts_name[4]].reshape([grid_lat,grid_lon])
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
	l3name  = output_prefix + '.A{:04d}{:02d}'.format(year,month)
	subname = '_baseline_monthly_v8.h5'
	ff=h5py.File(output_dir+l3name+subname,'w')

	PC=ff.create_dataset('lat_bnd',data=map_lat)
	PC.attrs['units']='degrees'
	PC.attrs['long_name']='Latitude_boundaries'    

	PC=ff.create_dataset('lon_bnd',data=map_lon)
	PC.attrs['units']='degrees'
	PC.attrs['long_name']='Longitude_boundaries'    

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
				addGridEntry(ff,new_name,unit_list[cnt],longname_list[cnt],fillvalue_list[cnt],scale_list[cnt],offst_list[cnt],grid_data[key])
				cnt += 1
	
	ff.close()

	print(l3name+subname+' Saved!')
	#---------------------------COMPLETED------------------------------------------------------

