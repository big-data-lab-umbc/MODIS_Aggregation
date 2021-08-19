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
			Change the definition of day to match with teh current C6 MYD08 product.
			Fixed the attributes in pixel count, histogram count and joint histogram count.
"""

import os 
import sys
import h5py
import timeit
import random
import calendar
import numpy as np
import pandas as pd 
from netCDF4 import Dataset
from collections import OrderedDict
from datetime import date, datetime
from dateutil.rrule import rrule, DAILY, MONTHLY

#This function is used for unbit the cloud mask from MYD06 product	
def bits_stripping(bit_start,bit_count,value):
    bitmask=pow(2,bit_start+bit_count)-1
    return np.right_shift(np.bitwise_and(value,bitmask),bit_start)

def read_filelist(MYD06_dir,MYD06_prefix,MYD03_dir,MYD03_prefix,year,day_in_year,time,fileformat):
	filename1,filename2 = [],[]

	for i in range(len(time)):

		if i >= 24: 
			yc ='%04i' % year[1]
			dc ='%03i' % day_in_year[1]
			tc ='%02i' % time[i]
		else:
			yc ='%04i' % year[0]
			dc ='%03i' % day_in_year[0]
			tc ='%02i' % time[i]

		# Read the filelist in the specific directory
		str_06 = os.popen("ls "+ MYD06_dir + MYD06_prefix + yc + dc + "."+ tc +"*."+fileformat).read()
		str_03 = os.popen("ls "+ MYD03_dir + MYD03_prefix + yc + dc + "."+ tc +"*."+fileformat).read()
		
		if len(str_06) == 0: 
			print("## ERROR!!! ")
			print("## The pool lacks of data in the defined time period becase of the new Definition of Day.")
			print("## Lack of data from Year: {}; Day in year: {}; Hours: 00:00 ~ 03:00 ".format(yc,dc))
			print("## Please either request a shorter period of the aggregation or add more data into the pool.")
			sys.exit()

		fname1 = np.array(str_06.split("\n"))
		fname1 = np.delete(fname1,len(fname1)-1)
		fname2 = np.array(str_03.split("\n"))
		fname2 = np.delete(fname2,len(fname2)-1)

		#print(fname1,fname2)

		if len(fname1) != len(fname2): 
			print("## ERROR!!! ")
			print("## The MYD06 is not consistent with MYD03 within a time (year-day-hour): {}-{}-{}".format(yc,dc,tc))
			print("## Either the MYD06 or MYD03 is lack of data (or has invalid data).")
			print("## Please check the Data_input_path in the data_path.csv, and make sure the input path has consistent data.")
			sys.exit()

		filename1 = np.append(filename1,fname1)
		filename2 = np.append(filename2,fname2)

	return filename1,filename2

def readEntry(key,ncf,spl_num) :
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
	
	#Sampling the variable if the data is 1km resolution 
	data_dim = len(str(rdval.shape[0]))
	if data_dim == 4: #1km resolution has 2030 row pixels, which is 4 digits of data_dim
		rdval = rdval[:,:1350] #Drop the last 4 pixels of the cross scan track if it is 1km product
		rdval = rdval[2::spl_num,3::spl_num]

	return rdval,data_dim,lonam,unit,fillvalue,scale,offst

def  read_MODIS(varnames,fname1,fname2, spl_num): 
	# Store the data from variables after reading MODIS files
	data={}
	
	# Read the User-defined variables from MYD06 product
	ncfile=Dataset(fname1,'r')
		
	# Read the User-defined input variables from MYD06 product
	for key in varnames:
		if key == 'cloud_fraction_CM': 
			continue #Ignoreing cloud_fraction from the input file
		else:
			data[key],data_dim,lonam,unit,fill,scale,offst = readEntry(key,ncfile,spl_num)
			data[key] = (data[key] - offst) / scale 
			data[key] = (data[key] - offst) * scale 
	
	# Read the Cloud Mask from MYD06 product
	if data_dim == 4:
		d06_CM = np.array(ncfile.variables['Cloud_Mask_1km'][:,:,0])
		attr_CM = ncfile.variables['Cloud_Mask_1km']._FillValue
		d06_CM[np.where(d06_CM == attr_CM)] = np.nan	
		d06_CM = d06_CM[:,:1350] #Drop the last 4 pixels of the cross scan track if it is 1km product
		d06_CM = d06_CM[2::spl_num,2::spl_num]
		if len(np.where(np.isnan(d06_CM) == 1)[0]) > 0: print("NaN cloud mask found!")
		CM1km  = np.array(d06_CM,dtype='byte') #[~np.isnan(d06_CM)]
	else:
		d06_CM = np.array(ncfile.variables['Cloud_Mask_5km'][:,:,0])
		attr_CM = ncfile.variables['Cloud_Mask_5km']._FillValue
		d06_CM[np.where(d06_CM == attr_CM)] = np.nan	
		if len(np.where(np.isnan(d06_CM) == 1)[0]) > 0: print("NaN cloud mask found!")
		CM1km  = np.array(d06_CM,dtype='byte')
		
	cloud_mask_flag = bits_stripping(1,2,CM1km)
	cloud_stat_flag = bits_stripping(0,1,CM1km)
	cloud_mask_flag[np.where(cloud_stat_flag == 0)] = -1 #exclude the non-determined pixels by the status flag of the first byte
	data['CM'] = cloud_mask_flag.astype(np.float)
	
	# This is for intepretation of cloud top properties QA from "Quality_Assurance_5km", which is not used in the current version
	#for key in varnames:
	#	if key == 'Cloud_Top_Pressure':
	#		CTP_useful = bits_stripping(0,1,QM1km)
	#		CTP_confid = bits_stripping(1,3,QM1km)
	#		CTP_QM = CTP_confid.astype(np.float)
	#		data[key][np.where(CTP_QM == 0)] = np.nan
	#	if key == 'Cloud_Top_Temperature':
	#		CTT_useful = bits_stripping(4,1,QM1km)
	#		CTT_confid = bits_stripping(5,3,QM1km)
	#		CTT_QM = CTT_confid.astype(np.float)
	#		data[key][np.where(CTT_QM == 0)] = np.nan
	#		#print(key,data[key].shape)
	#	if key == 'CM':
	#		CM_useful = bits_stripping(8,1,QM1km)
	#		CM_confid = bits_stripping(9,3,QM1km)
	#		CM_QM = CM_confid.astype(np.float)
	#		data[key][np.where(CM_QM == 0)] = 4
	#	if key == 'Cloud_Effective_Emissivity':
	#		CEE_useful = bits_stripping(12,1,QM1km)
	#		CEE_confid = bits_stripping(13,3,QM1km)
	#		CEE_QM = CEE_confid.astype(np.float)
	#		data[key][np.where(CEE_QM == 0)] = np.nan
		
	ncfile.close()

	# Read the common variables (Latitude & Longitude) 
	#   1km resolution: from MYD03 product 
	#	5km resolution: from MYD06 product
	
	if data_dim == 4:
		ncfile=Dataset(fname2,'r')
		d03_lat  = np.array(ncfile.variables['Latitude'][:,:])
		d03_lon  = np.array(ncfile.variables['Longitude'][:,:])
		d03_lat  = d03_lat[:,:1350] #Drop the last 4 pixels of the cross scan track if it is 1km product
		d03_lon  = d03_lon[:,:1350] #Drop the last 4 pixels of the cross scan track if it is 1km product
		lat  = d03_lat[2::spl_num,2::spl_num] #The geolocation is the 3th pixel cross and the 3th pixel along the swath
		lon  = d03_lon[2::spl_num,2::spl_num] #The geolocation is the 3th pixel cross and the 3th pixel along the swath
		attr_lat = ncfile.variables['Latitude']._FillValue
		attr_lon = ncfile.variables['Longitude']._FillValue
	else: 
		ncfile=Dataset(fname1,'r')
		lat  = np.array(ncfile.variables['Latitude'][:,:])
		lon  = np.array(ncfile.variables['Longitude'][:,:])
		attr_lat = ncfile.variables['Latitude']._FillValue
		attr_lon = ncfile.variables['Longitude']._FillValue

	# If the variable is not 1km product, exit and tell the User to reset the variables.
	for key in varnames:
		if key == 'cloud_fraction_CM': continue #Ignoreing cloud_fraction from the input file	
		if data[key].shape[0] != lat.shape[0]:
			print("## ERROR!!! ")
			print("## The dimension of varibale '"+key+"' is not match with latitude & longitude.")
			print("## Input variables should have the same resolution (either 1km or 5km).")
			print("## Check your varibales in files: \
				   {} \
				   {}  ".format(fname1,fname2))
			print(key+" shape:",data[key].shape,"Lat shape:",lat.shape)
			sys.exit()

	ncfile.close()

	return lat,lon,data


def cal_stats(z,key,grid_data,min_val,max_val,tot_val,count,all_val,all_val_2d, \
			  sts_switch,sts_name,intervals_1d,intervals_2d,key_idx, histnames):
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
	
	# Check if all of data within a grid is NaN, if yes, the STD and Histogram will not counted
	if len(all_val[~np.isnan(all_val)]) != 0: 
		
		grid_data['GRID_Counts'][z] += 1

		#Standard Deviation 
		if sts_switch[4] == True:
			if key == 'cloud_fraction_CM':
				grid_data[key+'_'+sts_name[4]][z] += tot_val**2
			else:
				grid_data[key+'_'+sts_name[4]][z] += np.sum(all_val[~np.isnan(all_val)]**2) #np.nansum((all_val - tot_val/count)**2) / count #tot_val**2
			
		#1D Histogram 
		if sts_switch[5] == True:	
			bin_interval1 = np.fromstring(intervals_1d[key_idx], dtype=np.float, sep=',' )
			if all_val.size == 1: 
				all_val = np.array([all_val])
			else:
				hist_idx1 = np.histogram(all_val[~np.isnan(all_val)],bins=bin_interval1)[0]
				grid_data[key+'_'+sts_name[5]][z,:] += hist_idx1
			
		#2D Histogram 
		if sts_switch[6] == True:
			bin_interval1 = np.fromstring(intervals_1d[key_idx], dtype=np.float, sep=',' )
			bin_interval2 = np.fromstring(intervals_2d[key_idx], dtype=np.float, sep=',' )
			if all_val.size == 1:
				all_val = np.array([all_val])
				all_val_2d = np.array([all_val_2d])
			else:
				if all_val[~np.isnan(all_val)].shape[0] != all_val_2d[~np.isnan(all_val_2d)].shape[0]:
					grid_data[key+'_'+sts_name[6]+histnames[key_idx]][z,:,:] += 0
				else:
					hist_idx2 = np.histogram2d(all_val[~np.isnan(all_val)],all_val_2d[~np.isnan(all_val_2d)],bins=(bin_interval1,bin_interval2))[0]
					grid_data[key+'_'+sts_name[6]+histnames[key_idx]][z,:,:] += hist_idx2

		
	return grid_data

def run_modis_aggre(fname1,fname2,day_in_year,shift_hour,NTA_lats,NTA_lons,grid_lon,grid_lat,gap_x,gap_y,hdfs, \
					grid_data,sts_switch,varnames,intervals_1d,intervals_2d,var_idx, spl_num, sts_name, histnames):
	# This function is the data aggregation loops by number of files
	hdfs = np.array(hdfs)
	for j in hdfs: #range(2): range(143,145): 
		print("File Number: {} / {}".format(j,hdfs[-1]))
		
		# Retrieve the day and hour from the file name
		file_day  = np.int(fname1[j].split('.')[1][5:])
		file_hour = np.int(fname1[j].split('.')[2][:2])
		
		# Read Level-2 MODIS data
		lat,lon,data = read_MODIS(varnames,fname1[j],fname2[j],spl_num)
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
			if key == 'cloud_fraction_CM': 
				CF_key_idx = key_idx
				key_idx += 1
				continue #Ignoreing Cloud_Fraction from the input file	
			data[key] = data[key][res_idx].ravel()
			key_idx += 1
		
		# Artifact Adjustment of definition of day
		# Daytime in longitude from -90 to -180
		if (file_day == day_in_year[0]) & (file_hour < shift_hour): #& (np.nanmean(np.abs(lon)) > 145) :
			cutoff = np.where((lon <= -90) & (lon >= -180)) #(lon < 0)
			for key in varnames:
				if key == 'cloud_fraction_CM':
					CM[cutoff] =np.nan
					continue
				data[key][cutoff] = np.nan
			
		if (file_day == day_in_year[1]) & (file_hour < shift_hour): #& (np.nanmean(np.abs(lon)) > 145) :
			cutoff = np.where((lon >= 90) & (lon <= 180)) #(lon-180) > -180
			for key in varnames:
				if key == 'cloud_fraction_CM':
					CM[cutoff] =np.nan
					continue
				data[key][cutoff] = np.nan

		# Nighttime in longitude from 0 to 90
		if (file_day == day_in_year[0]) & (file_hour < shift_hour): #(np.nanmean(np.abs(lon)) < 35) :
			cutoff = np.where((lon >= 0) & (lon <= 90))
			for key in varnames:
				if key == 'cloud_fraction_CM':
					CM[cutoff] =np.nan
					continue
				data[key][cutoff] = np.nan

		if (file_day == day_in_year[1]) & (file_hour < shift_hour): #(np.nanmean(np.abs(lon)) < 35) :
			cutoff = np.where((lon <= 0) & (lon >= -90))
			for key in varnames:
				if key == 'cloud_fraction_CM':
					CM[cutoff] =np.nan
					continue
				data[key][cutoff] = np.nan

		# Locate the lat lon index into 3-Level grid box
		idx_lon = ((lon-NTA_lons[0])/gap_x).astype(int)
		idx_lat = ((lat-NTA_lats[0])/gap_y).astype(int) 

		latlon_index=(idx_lat*grid_lon)+idx_lon

		latlon_index_unique = np.unique(latlon_index)

		for i in np.arange(latlon_index_unique.size):
		#-----loop through all the grid boxes ocupied by this granule------#
			z=latlon_index_unique[i]
			if((z >= 0) & (z < (grid_lat*grid_lon))):

				# For cloud fraction
				TOT_pix = np.sum((CM[np.where(latlon_index == z)]>=0)&(CM[np.where(latlon_index == z)]<=3)).astype(float)
				CLD_pix = np.sum((CM[np.where(latlon_index == z)]>=0)&(CM[np.where(latlon_index == z)]<=1)).astype(float)
				
				Fraction = CLD_pix / TOT_pix
				
				if len(intervals_2d) != 1:	
					pixel_data_2d = data[varnames[var_idx[CF_key_idx]]]
					ave_val_2d = np.nansum(pixel_data_2d[np.where(latlon_index == z)]).astype(float) / TOT_pix
				else:
					ave_val_2d = 0

				# Calculate Statistics pamameters
				grid_data = cal_stats(z,"cloud_fraction_CM",grid_data, \
									  Fraction,Fraction,CLD_pix,TOT_pix,Fraction,ave_val_2d, \
									  sts_switch,sts_name,intervals_1d,intervals_2d,CF_key_idx, histnames)

				# For other variables
				key_idx = 0
				for key in varnames:
					if key == 'cloud_fraction_CM': #Ignoreing cloud_fraction from the input file	
						key_idx += 1
						continue 

					pixel_data = data[key]
					
					tot_val = np.nansum(pixel_data[np.where(latlon_index == z)]).astype(float)
					TOT_pix = np.where(~np.isnan(pixel_data[np.where(latlon_index == z)]))[0].shape[0]
					all_val = np.array(pixel_data[np.where(latlon_index == z)]).astype(float)
					max_val = np.nanmax(pixel_data[np.where(latlon_index == z)]).astype(float)
					min_val = np.nanmin(pixel_data[np.where(latlon_index == z)]).astype(float)

					if len(intervals_2d) != 1:
						pixel_data_2d = data[varnames[var_idx[key_idx]]]
						all_val_2d = np.array(pixel_data_2d[np.where(latlon_index == z)]).astype(float)
					else:
						all_val_2d = 0

					# Calculate Statistics pamameters
					grid_data = cal_stats(z,key,grid_data, \
										  min_val,max_val,tot_val,TOT_pix,all_val,all_val_2d, \
										  sts_switch,sts_name,intervals_1d,intervals_2d,key_idx, histnames)

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
		scale_factor = 1.0
		add_offset = 0.0
	elif (('Maximum' in name) == True) | (('Minimum' in name) == True):
		tmp_data = data/scale_factor + add_offset
		tmp_data[np.where(np.isinf(tmp_data) == 1)]=fillvalue
		original_data = tmp_data.astype(np.int)
	else: 
		tmp_data = data/scale_factor + add_offset
		tmp_data[np.where(np.isnan(tmp_data) == 1)]=fillvalue
		original_data = tmp_data.astype(np.int)

	# This function is to adjust the scaling factor of cloud fraction if is computed by the 5-km cloud fraction from MOD06/MYD06
	if (('Cloud_Fraction' in name) == True):
		scale_factor = 0.0001
		add_offset = 0.0
		fillvalue = -9999
		tmp_data = data/scale_factor + add_offset
		tmp_data[np.where(np.isnan(tmp_data) == 1)]=fillvalue
		original_data = tmp_data

	PCentry=f.create_dataset(name,data=original_data)
	PCentry.dims[0].label='lat_bnd'
	PCentry.dims[1].label='lon_bnd'
	PCentry.attrs['units']       = np.str(units)
	PCentry.attrs["long_name"]   = np.str(long_name)
	PCentry.attrs['_FillValue']  = fillvalue
	PCentry.attrs['scale_factor']= scale_factor
	PCentry.attrs['add_offset']  = add_offset

def read_user_inputs():

	#-------------STEP 1: Get all the inputs from User --------
	# checking user input
	if (len(sys.argv) != 16) & (len(sys.argv) != 17):
		print("## ERROR!!!")
		print("## Wrong user input")
		print("## Please see the example below:")
		print("## python MODIS_Aggregation_Serial.py \
						 <Data Path> <Start Date: yyyy-mm-dd> <End Date: yyyy-mm-dd> \
						 <Polygon boundaries> <Lat & Lon Grid Size > \
						 <Sampling number larger than 0> \
						 <1/0> <1/0> <1/0> \
						 <1/0> <1/0> <1/0> \
						 <1/0> <Variable Imput File> <JHist Variable Imput File>")
		
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
		
	#-------------STEP 2: Set up the specific directory --------
	data_path_file = np.array(pd.read_csv(sys.argv[1], header=0, delim_whitespace=True))
	MYD06_dir    = data_path_file[0,0] #'/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD06_L2/'
	MYD06_prefix = data_path_file[0,1] #'MYD06_L2.A'
	MYD03_dir    = data_path_file[1,0] #'/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD03/'
	MYD03_prefix = data_path_file[1,1] #'MYD03.A'
	fileformat = 'hdf'

	output_path_file = np.array(pd.read_csv(sys.argv[1], header=3, delim_whitespace=True))
	output_dir = output_path_file[0,0]
	output_prefix = output_path_file[0,1]

	#-------------STEP 3: Set up spactial and temporal resolution & variable names----------
	NTA_lats = [poly[0],poly[1]] 
	NTA_lons = [poly[2],poly[3]] 
	
	gap_x, gap_y = grid[1],grid[0] 

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

	#--------------STEP 4: Create arrays for level-3 statistics data-------------------------
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

	#--------------STEP 5: Read the filename list for 3-hours time shifting due to the new definition of day-------------------
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
				
	filenum = np.arange(len(fname1))
	
	#--------------STEP 6:  Check if output hdf file name already exist
	l3name  = output_prefix + '.A{:04d}{:03d}-{:04d}{:03d}.h5'.format(year[0],day_in_year[0],year[1],day_in_year[1])
	hdf_exist = np.int(os.popen('[ -f {} ] && echo "1" || echo "2"'.format(output_dir+l3name)).read()[:-1])
	if hdf_exist == 1: 
		print("## ERROR!!!")
		print("## The output HDF file name exist")
		print("## Please rename the output HDF!")
		sys.exit()
	
	#--------------STEP 7: Read Attributes of each variables----------------------------------
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
		if key == 'cloud_fraction_CM': 
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

	# Add the long name of the cloud freaction statistic at the first row
	CM_unit         = 'none'
	CM_longname     = 'Cloud Fraction from Cloud Mask Flag (cloudy & prob cloudy)'
	CM_fillvalue    = -9999
	CM_scale_factor = 0.0001
	CM_add_offset   = 0.0
	unit_list      = np.insert(unit_list,      name_idx, CM_unit)
	scale_list     = np.insert(scale_list,     name_idx, CM_scale_factor)
	offst_list     = np.insert(offst_list,     name_idx, CM_add_offset)
	longname_list  = np.insert(longname_list,  name_idx, CM_longname)
	fillvalue_list = np.insert(fillvalue_list, name_idx, CM_fillvalue)

	ncfile.close()

	return fname1,fname2,day_in_year,shift_hour,NTA_lats,NTA_lons,map_lon,map_lat,grid_lon,grid_lat,gap_x,gap_y,filenum, \
		   grid_data,sts_switch,varnames,intervals_1d,intervals_2d,bin_num1,bin_num2,var_idx,spl_num,sts_name,histnames, \
		   output_dir,l3name,unit_list,scale_list,offst_list,longname_list,fillvalue_list

if __name__ =='__main__':
# This is the main program for using concurrent to speed up the whole process
	
	#--------------STEP 1: Read User Inputs and Initial Paramters for Aggregation--------------------
	fname1,fname2,day_in_year,shift_hour,NTA_lats,NTA_lons,map_lon,map_lat,grid_lon,grid_lat,gap_x,gap_y,filenum, \
	grid_data,sts_switch,varnames,intervals_1d,intervals_2d,bin_num1,bin_num2,var_idx,spl_num,sts_name,histnames, \
	output_dir,l3name,unit_list,scale_list,offst_list,longname_list,fillvalue_list = read_user_inputs()

	#--------------STEP 2: Start Internal Aggregation------------------------------------------------

	# Start counting operation time
	start_time = timeit.default_timer() 

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
				if sts_idx[i] >= 5:
					addGridEntry(ff,new_name,unit_list[cnt],longname_list[cnt],fillvalue_list[cnt],scale_list[cnt],offst_list[cnt],grid_data[key],intervals_1d[cnt],intervals_2d[cnt])
				else:
					addGridEntry(ff,new_name,unit_list[cnt],longname_list[cnt],fillvalue_list[cnt],scale_list[cnt],offst_list[cnt],grid_data[key],intervals_1d[0],intervals_2d[0])
				cnt += 1
	
	ff.close()

	print(l3name+' Saved!')
	#---------------------------COMPLETED------------------------------------------------------

