#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Main Program: RUN MODIS AGGREGATION IN PARALLEL ON FILE LEVEL

Created on 2019

@author: Jianyu Zheng
		 jzheng3@umbc.edu
		 Atmospheric Physics
		 University of Maryland, Baltimore County
"""

import os 
import sys
import timeit
import random
import numpy as np
from mpi4py import MPI
from netCDF4 import Dataset

def read_filelist(loc_dir,prefix,unie,fileformat):
	# Read the filelist in the specific directory
	str = os.popen("ls "+ loc_dir + prefix + unie + "*."+fileformat).read()
	fname = np.array(str.split("\n"))
	fname = np.delete(fname,len(fname)-1)

	return fname

def read_MODIS(fname1,fname2,verbose=False): # READ THE HDF FILE
	
	# Read the cloud mask from MYD06_L2 product')
	ncfile=Dataset(fname1,'r')
	CM1km = np.array(ncfile.variables['Cloud_Mask_1km'])
	CM   = (np.array(CM1km[:,:,0],dtype='byte') & 0b00000110) >>1
	ncfile.close()

	# Read the geolocation data from MYD03 product')
	ncfile=Dataset(fname2,'r')
	lat  = np.array(ncfile.variables['Latitude'])
	lon  = np.array(ncfile.variables['Longitude'])
	attr_lat = ncfile.variables['Latitude']._FillValue
	attr_lon = ncfile.variables['Longitude']._FillValue

	#Use _FillValue to remove fill data in lat & lon
	lat[np.where(lat == attr_lat)] = 0.0
	lon[np.where(lat == attr_lat)] = 0.0
	CM [np.where(lat == attr_lat)] = 0.5 #which will not be identified by lines 80-83 

	lat[np.where(lon == attr_lon)] = 0.0
	lon[np.where(lon == attr_lon)] = 0.0
	CM [np.where(lon == attr_lon)] = 0.5 #which will not be identified by lines 80-83
	ncfile.close()

	return lat,lon,CM

def run_modis_aggre(TOT_pix,CLD_pix,hdfs):
	# This function is the data aggregation loops by number of files
	hdfs = np.array(hdfs)

	if hdfs[0] < 31: 

		dc ='%03i' % hdfs[0]

		fname1 = read_filelist(MYD06_dir,MYD06_prefix,dc,fileformat)
		fname2 = read_filelist(MYD03_dir,MYD03_prefix,dc,fileformat)

		for j in np.arange(len(fname1)):
	
			#print("File Number: {} / {}".format(j,len(fname1)))
			
			# Read Level-2 MODIS data
			lat,lon,CM = read_MODIS(fname1[j],fname2[j])
			#print('CM:',CM)
			# Ravel the 2-D data to 1-D array
			lat = lat.ravel()
			lon = lon.ravel()
			CM  = CM.ravel()
				
			lat = lat.astype(int)
			lon = lon.astype(int)
			latlon_index=((lat+90)*grid_lon)+(lon+180)
		
			latlon_index_unique = np.unique(latlon_index)
			for i in np.arange(latlon_index_unique.size):
			#-----loop through all the grid boxes ocupied by this granule------#
				z=latlon_index_unique[i]
				if( z >= 0):
					TOT_pix[z] = np.sum(CM[np.where(latlon_index == z)]>=0)
					CLD_pix[z] = np.sum(CM[np.where(latlon_index == z)]<=1)
	else:
		TOT_pix[:] = 0
		CLD_pix[:] = 0

	#print('File name:',CLD_pix[40:50],TOT_pix[40:50])
	#print(TOT_pix.shape,latlon_index_unique,latlon_index_unique.shape)
	#fout1 = open('test.dat','w')
	#for i in range(len(TOT_pix)):
	#	fout1.write('{:<8.3f}\n'.format(TOT_pix[i]))
	#fout1.close()
	#print(TOT_pix)
	
	return (TOT_pix,CLD_pix)	

##Run the code in series for days in one month
#def run_modis(day_num):
#	days = np.linspace(1,day_num,day_num)
#	for day in days:
#		day_fraction=run_modis_aggre(day)
#	return day_fraction

if __name__ =='__main__':
# This is the main program for using concurrent to speed up the whole process
	
	#-------------STEP 1: Set up the specific directory --------
	MYD06_dir= '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD06_L2/'
	MYD06_prefix = 'MYD06_L2.A2008'
	MYD03_dir= '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD03/'
	MYD03_prefix = 'MYD03.A2008'
	fileformat = 'hdf'
	
	#-------------STEP 2: Set up spactial and temporal resolution----------
	NTA_lats = [-90,90]   
	NTA_lons = [-180,180] 
	
	gap_x, gap_y = 1, 1
	
	grid_lon=np.int((NTA_lons[-1]-NTA_lons[0])/gap_x)
	grid_lat=np.int((NTA_lats[-1]-NTA_lats[0])/gap_y)
	
	TOT_pix     = np.zeros(grid_lat*grid_lon)
	CLD_pix     = np.zeros(grid_lat*grid_lon)
	sum_tot_pix = np.zeros(grid_lat*grid_lon)
	sum_cld_pix = np.zeros(grid_lat*grid_lon)

	# Initiate MPI 
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	random.seed(rank)
	#print(size)

	# Initiate the number of day and total cloud fraction
	day_num   = np.linspace(1,32,32,dtype=np.int)

	# Distribute the day's loops into ppns
	tasks = np.array(np.split(day_num,size))
	hdfs = tasks[rank]
	
	print("process {} aggregating days {} / {}...".format(rank, hdfs[0],day_num[-1]))
	
	start_time = timeit.default_timer() #time.perf_counter() # Start counting operation time
	results = np.asarray(run_modis_aggre(TOT_pix,CLD_pix,hdfs))
		
	if rank == 0:
		sum_tot_pix = sum_tot_pix + results[0,:]
		sum_cld_pix = sum_cld_pix + results[1,:]
		for i in range(1,size):
			results = comm.recv(source=i, tag=0)
			#recv_req = comm.Irecv(results,source=i, tag=0)
			#recv_req.wait()
			
			sum_tot_pix = sum_tot_pix + results[0,:] #TOT_pix
			sum_cld_pix = sum_cld_pix + results[1,:] #CLD_pix

		# Compute the monthly average cloud fraction 
		ave_fra_grid = (sum_cld_pix / sum_tot_pix).reshape([grid_lat,grid_lon])
		end_time = timeit.default_timer() #time.perf_counter() # End counting operation time
		print('ave_fra_grid:')
		print( ave_fra_grid  )

		print ("Operation Time in {:7.2f} seconds".format(end_time - start_time))

		# Create file to store the result 
		np.savetxt("cld_fraction_mpi_days.dat", ave_fra_grid, fmt="%10.4f")
		
	else:
		print("Process {} finished".format(rank))
		comm.send(results, dest=0, tag=0)
		#send_req = comm.Isend(results, dest=0, tag=0)
		#send_req.wait()
