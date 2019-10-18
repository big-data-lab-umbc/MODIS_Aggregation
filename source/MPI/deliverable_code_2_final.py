#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Main Program: Run MODIS AGGREGATION IN PARALLEL

Created on 2019

@author: Jianyu Zheng
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

def run_modis_aggre(fname1,fname2,NTA_lats,NTA_lons,grid_lon,gap_x,gap_y,TOT_pix,CLD_pix,hdfs):
	# This function is the data aggregation loops by number of files
	hdfs = np.array(hdfs)
	for j in hdfs:
		print("File Number: {} / {}".format(j,hdfs[-1]))
	
		# Read Level-2 MODIS data
		lat,lon,CM = read_MODIS(fname1[j],fname2[j])
		#print(lat.shape,lon.shape,CM.shape)

		# Restrain lat & lon & variables in the required region 
		res_idx = np.where((lat > NTA_lats[0]) & (lat < NTA_lats[1]) & (lon > NTA_lons[0]) & (lon < NTA_lons[1]))
		#print(res_idx)
		CM  = CM [res_idx]
		lat = lat[res_idx]
		lon = lon[res_idx]

		# Ravel the 2-D data to 1-D array
		lat = lat.ravel()
		lon = lon.ravel()
		CM  = CM.ravel()
		
		#if len(CM) >1: print(lat[len(CM)-1]) #18496, 18496, 16199, 785409, 28.731682, 22.853891)

		# Locate the lat lon index into 3-Level frid box
		idx_lon = ((lon-NTA_lons[0])/gap_x).astype(int)
		idx_lat = ((lat-NTA_lats[0])/gap_y).astype(int)

		latlon_index=(idx_lat*grid_lon)+idx_lon

		latlon_index_unique = np.unique(latlon_index)

		for i in np.arange(latlon_index_unique.size):
		#-----loop through all the grid boxes ocupied by this granule------#
			z=latlon_index_unique[i]
			if((z >= 0) & (z < len(TOT_pix))):
				TOT_pix[z] = np.sum(CM[np.where(latlon_index == z)]>=0)
				CLD_pix[z] = np.sum(CM[np.where(latlon_index == z)]<=1)

	return (TOT_pix,CLD_pix)	

if __name__ =='__main__':
# This is the main program for using concurrent to speed up the whole process
	
	#-------------STEP 1: Set up the specific directory --------
	MYD06_dir= '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD06_L2/'
	MYD06_prefix = 'MYD06_L2.A2008'
	MYD03_dir= '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD03/'
	MYD03_prefix = 'MYD03.A2008'
	fileformat = 'hdf'
	
	#-------------STEP 2: Set up spactial and temporal resolution----------
	NTA_lats = [-90,90]   #[  0,40] #[-90,90]   #[-30,30]    
	NTA_lons = [-180,180] #[-40,60] #[-180,180] #[-60,60]  
	
	gap_x, gap_y = 1,1 #0.625,0.5

	if ((NTA_lons[-1]-NTA_lons[0])%gap_x != 0) | ((NTA_lats[-1]-NTA_lats[0])%gap_y != 0): 
		print("Grid size should be dividable by the dimension of the selected region. Please try again!")
		sys.exit()

	map_lon = np.arange(NTA_lons[0],NTA_lons[1],gap_x)
	map_lat = np.arange(NTA_lats[0],NTA_lats[1],gap_y)
	Lon,Lat = np.meshgrid(map_lon,map_lat)
	
	grid_lon=np.int((NTA_lons[-1]-NTA_lons[0])/gap_x)
	grid_lat=np.int((NTA_lats[-1]-NTA_lats[0])/gap_y)

	#print(grid_lon,grid_lat,grid_lat*grid_lon)
	
	TOT_pix     = np.zeros(grid_lat*grid_lon)
	CLD_pix     = np.zeros(grid_lat*grid_lon)
	sum_tot_pix = np.zeros(grid_lat*grid_lon)
	sum_cld_pix = np.zeros(grid_lat*grid_lon)

	fname1,fname2 = [],[]

	# Read all files in a month (in this case: January)
	days = np.arange(1,2,dtype=np.int)
	for day in days:
		dc ='%03i' % day
		fname_tmp1 = read_filelist(MYD06_dir,MYD06_prefix,dc,fileformat)
		fname_tmp2 = read_filelist(MYD03_dir,MYD03_prefix,dc,fileformat)
		fname1 = np.append(fname1,fname_tmp1)
		fname2 = np.append(fname2,fname_tmp2)

	# Initiate MPI 
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	random.seed(rank)

	# Distribute the number of files into ppns for MPI
	remain   = size-len(fname1)%size
	ppn_file = (len(fname1)+remain)/size 

	if ppn_file > remain: 
		# Distribute the day's loops into MPI ppns
		files = np.arange(len(fname1)+remain)
		tasks = np.array(np.split(files,size))
		hdfs = tasks[rank]
	
		if rank == (size-1): 
			hdfs = np.delete(hdfs, np.arange(len(hdfs)-remain,len(hdfs)))
	else:
		# Distribute the day's loops into MPI ppns
		files = np.arange(len(fname1)-len(fname1)%size)
		tasks = np.array(np.split(files,size))
		hdfs = tasks[rank]
	
		if rank == (size-1): 
			hdfs = np.append(hdfs, np.arange(len(files),len(files)+len(fname1)%size))

	print("process {} aggregating files from {} to {}...".format(rank, hdfs[0],hdfs[-1]))
	
	# Start counting operation time
	start_time = timeit.default_timer() 

	results = np.asarray(run_modis_aggre(fname1,fname2,NTA_lats,NTA_lons,grid_lon,gap_x,gap_y,TOT_pix,CLD_pix,hdfs))
		
	if rank == 0:
		TOT_pix = TOT_pix + results[0,:]
		CLD_pix = CLD_pix + results[1,:]
		for i in range(1,size):
			recv_req = comm.Irecv(results,source=i, tag=0)
			recv_req.wait()
			
			TOT_pix = TOT_pix + results[0,:] #TOT_pix
			CLD_pix = CLD_pix + results[1,:] #CLD_pix

		# Compute the monthly average cloud fraction 
		ave_fra_grid = (CLD_pix / TOT_pix).reshape([grid_lat,grid_lon])

		end_time = timeit.default_timer()

		print('ave_fra_grid:')
		print( ave_fra_grid  )

		print ("Operation Time in {:7.2f} seconds".format(end_time - start_time))
		
		# Create file to store the result 
		np.savetxt("test_cloud_fraction.dat"  , ave_fra_grid, fmt="%10.4f")
		np.savetxt("test_geolocation_lat.dat" , Lat, fmt="%10.4f")
		np.savetxt("test_geolocation_lon.dat" , Lon, fmt="%10.4f")
		
	else:
		print("Process {} finished".format(rank))
		send_req = comm.Isend(results, dest=0, tag=0)
		send_req.wait()