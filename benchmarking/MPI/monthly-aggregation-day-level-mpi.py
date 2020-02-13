#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Main Program: Run MODIS AGGREGATION IN PARALLEL BY MPI

Created on 2019

@author: Jianyu Zheng
"""

import os 
import sys
import h5py
import glob
import timeit
import random
import numpy as np
import xarray as xr
from mpi4py import MPI
from netCDF4 import Dataset
import matplotlib.pyplot as plt

def read_MODIS(M06_files,M03_files,verbose=False): # READ THE HDF FILE
	
	# Read the cloud mask from MYD06_L2 product')
	ncfile=Dataset(M06_files,'r')
	d06_CM = ncfile.variables['Cloud_Mask_1km'][:,:,0]
	CM1km  = d06_CM[::3,::3]
	CM     = (np.array(CM1km,dtype='byte') & 0b00000110) >>1
	ncfile.close()

	# Read the geolocation data from MYD03 product')
	ncfile=Dataset(M03_files,'r')
	d03_lat  = ncfile.variables['Latitude'][:,:]
	d03_lon  = ncfile.variables['Longitude'][:,:]
	lat  = d03_lat[::3,::3]
	lon  = d03_lon[::3,::3]
	attr_lat = ncfile.variables['Latitude']._FillValue
	attr_lon = ncfile.variables['Longitude']._FillValue

	#Use _FillValue to remove fill data in lat & lon
	#lat[np.where(lat == attr_lat)] = 0.0
	#lon[np.where(lat == attr_lat)] = 0.0
	#CM [np.where(lat == attr_lat)] = 0.5 #which will not be identified by lines 80-83 
#
	#lat[np.where(lon == attr_lon)] = 0.0
	#lon[np.where(lon == attr_lon)] = 0.0
	#CM [np.where(lon == attr_lon)] = 0.5 #which will not be identified by lines 80-83
	ncfile.close()

	return lat,lon,CM

def run_modis_aggre(dayloop):
	# This function is the data aggregation loops by number of files
	dayloop = np.array(dayloop)+1

	for day in dayloop:

		if day > 31: break 
		
		dc ='%03i' % day
		M03_files = sorted(glob.glob(MYD03_dir + "MYD03.A2008" + dc + "*"))
    	M06_files = sorted(glob.glob(MYD06_dir + "MYD06_L2.A2008" + dc + "*"))

		for j in range(len(M06_files)):
			print("File Number: {} / {} in day {}".format(j,len(M06_files),day))
			
			# Read Level-2 MODIS data
			lat,lon,CM = read_MODIS(M06_files[j],M03_files[j])
			#print(lat.shape,lon.shape,CM.shape)
	
			# Ravel the 2-D data to 1-D array
			lat = (lat.ravel()+ 89.5).astype(int)
			lon = (lon.ravel()+ 179.5).astype(int)
			lat = np.where(lat > -1, lat, 0)
			lon = np.where(lon > -1, lon, 0)
			
			# increment total_pix by 1 for the grid for each value in (lat, lon).
			for i, j in zip(lat, lon):
					total_pix[i,j] += 1
			
			# covert ds06_decoded from 2D to 1D, check whether each element is less than or equal to 0, return a tuple whose first element is an 1D arrays of indices of ds06_decoded's elements whose value is less than or equal to 0.  
			index = np.nonzero(CM.ravel() == 0)
			# get its lat and lon for each cloud pixel.
			# we can use this approach because the internal structure (677, 452) is the same for both MYD03 and MYD06.
			cloud_lon = [lon[i] for i in index[0]]
			cloud_lat = [lat[i] for i in index[0]]
			 # increment cloud_pix by 1 for the grid for each value in (cloud_lat, cloud_lon).
			for x, y in zip(cloud_lat, cloud_lon):
				cloud_pix[x,y] += 1

	return (total_pix,cloud_pix)

def save_output(cf):
	cf1 = xr.DataArray(cf)
	cf1.to_netcdf("monthlyCloudFraction-day-level-parallelization.nc")
	plt.figure(figsize=(14, 7))
	plt.contourf(range(-180, 180), range(-90, 90), cf, 100, cmap="jet")
	plt.xlabel("Longitude", fontsize=14)
	plt.ylabel("Latitude", fontsize=14)
	plt.title("Level 3 Cloud Fraction Aggregation for January 2008", fontsize=16)
	plt.colorbar()
	plt.savefig("monthlyCloudFraction-day-level-parallelization.png")


if __name__ =='__main__':
# This is the main program for using concurrent to speed up the whole process

	# Start counting operation time
	start_time = timeit.default_timer()
			
	#-------------STEP 1: Define Files Path--------
	MYD06_dir= '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD06_L2/'
	MYD03_dir= '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD03/'
	
	#-------------STEP 2: Set up spactial and temporal resolution----------
	total_pix = np.zeros((180, 360))
	cloud_pix = np.zeros((180, 360))

	# Initiate MPI 
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	random.seed(rank)

	# Distribute the number of files into ppns for MPI
	day_num   = np.linspace(1,31,31,dtype=np.int)

	remain   = size-len(day_num)%size
	ppn_file = (len(day_num)+remain)/size 

	if len(day_num) <= size: 
		files = np.arange(len(day_num)+remain)
		tasks = np.array(np.split(files,size))
		dayloop = tasks[rank]

	elif ppn_file >= remain: 
		# Distribute the day's loops into MPI ppns
		files = np.arange(len(day_num)+remain)
		tasks = np.array(np.split(files,size))
		dayloop = tasks[rank]
	
		if rank == (size-1): 
			dayloop = np.delete(dayloop, np.arange(len(dayloop)-remain,len(dayloop)))
	else:
		# Distribute the day's loops into MPI ppns
		files = np.arange(len(day_num)-len(day_num)%size)
		tasks = np.array(np.split(files,size))
		dayloop = tasks[rank]
	
		if rank == (size-1): 
			dayloop = np.append(dayloop, np.arange(len(files),len(files)+len(day_num)%size))

	print("process {} aggregating days from {} to {}...".format(rank, dayloop[0],dayloop[-1]))
	
	# Start counting operation time
	# start_time = timeit.default_timer() 

	results = np.asarray(run_modis_aggre(dayloop))
		
	if rank == 0:
		total_pix += results[0,:]
		cloud_pix += results[1,:]
		
		for i in range(1,size):
			recv_req = comm.Irecv(results,source=i, tag=0)
			recv_req.wait()
		
			total_pix += results[0,:]
			cloud_pix += results[1,:]

		# Compute the mean cloud fraction & Statistics (Include Min & Max & Standard deviation)
		Mean_Fraction = (total_pix / cloud_pix)

		end_time = timeit.default_timer()
		
		# Create HDF5 file to store the result 
		save_output(Mean_Fraction)

		# end_time = timeit.default_timer()

		print('Mean_Fraction:')
		print( Mean_Fraction  )

		print ("Operation Time in {:7.2f} seconds".format(end_time - start_time))
		

	else:
		print("Process {} finished".format(rank))
		send_req = comm.Isend(results, dest=0, tag=0)
		send_req.wait()
