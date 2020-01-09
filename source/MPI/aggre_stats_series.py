#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Main Program: Run MODIS AGGREGATION IN PARALLEL

Created on 2019

@author: Jianyu Zheng
"""

#Updates: Add statistics for flexible variables

import os 
import sys
import h5py
import timeit
import random
import numpy as np
from mpi4py import MPI
from netCDF4 import Dataset

def read_filelist(loc_dir,prefix,yr,day,fileformat):
	# Read the filelist in the specific directory
	str = os.popen("ls "+ loc_dir + prefix + yr + day + "*."+fileformat).read()
	fname = np.array(str.split("\n"))
	fname = np.delete(fname,len(fname)-1)

	return fname

def readEntry(key,ncf):
	# Read the MODIS variables based on User's name list
    rdval=np.array(ncf.variables[key])
    scale=ncf.variables[key].scale_factor
    offst=ncf.variables[key].add_offset
    return (rdval+offst)*scale

def read_MODIS(varnames,fname1,fname2): 
	# Store the data from variables after reading MODIS files
	data={}

	# Read the common variables (Cloud Mask) from MYD06 product
	ncfile=Dataset(fname1,'r')
	CM1km = readEntry('Cloud_Mask_1km',ncfile)
	CM1km = np.array(ncfile.variables['Cloud_Mask_1km'])
	data['CM'] = (np.array(CM1km[:,:,0],dtype='byte') & 0b00000110) >>1
	for key in varnames:
		print(key)
		data[key] = readEntry(key,ncfile)
	ncfile.close()

	# Read the common variables (Latitude & Longitude) from MYD03 product
	ncfile=Dataset(fname2,'r')
	lat  = np.array(ncfile.variables['Latitude'])
	lon  = np.array(ncfile.variables['Longitude'])
	attr_lat = ncfile.variables['Latitude']._FillValue
	attr_lon = ncfile.variables['Longitude']._FillValue

	#Use _FillValue to remove fill data in lat & lon
	lat[np.where(lat == attr_lat)] = np.nan
	lon[np.where(lat == attr_lat)] = np.nan
	data['CM'][np.where(lat == attr_lat)] = np.nan #which will not be identified by lines 80-83 

	lat[np.where(lon == attr_lon)] = np.nan
	lon[np.where(lon == attr_lon)] = np.nan
	data['CM'][np.where(lon == attr_lon)] = np.nan #which will not be identified by lines 80-83
	ncfile.close()

	return lat,lon,data

def run_modis_aggre(fname1,fname2,NTA_lats,NTA_lons,grid_lon,gap_x,gap_y,hdfs, \
					grid_data,sts_switch,varnames,bin_num,lobnd1,upbnd1,lobnd2,upbnd2):
	# This function is the data aggregation loops by number of files
	hdfs = np.array(hdfs)
	for j in hdfs:
		print("File Number: {} / {}".format(j,hdfs[-1]))
	
		# Read Level-2 MODIS data
		lat,lon,data = read_MODIS(varnames,fname1[j],fname2[j])
		CM = data['CM']

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
		
		# Locate the lat lon index into 3-Level frid box
		idx_lon = ((lon-NTA_lons[0])/gap_x).astype(int)
		idx_lat = ((lat-NTA_lats[0])/gap_y).astype(int)

		latlon_index=(idx_lat*grid_lon)+idx_lon

		latlon_index_unique = np.unique(latlon_index)

		for i in np.arange(latlon_index_unique.size):
		#-----loop through all the grid boxes ocupied by this granule------#
			z=latlon_index_unique[i]
			if((z >= 0) & (z < len(Count))):
				TOT_pix = np.sum(CM[np.where(latlon_index == z)]>=0).astype(float)
				CLD_pix = np.sum(CM[np.where(latlon_index == z)]<=1).astype(float)

				Fraction = (CLD_pix / TOT_pix)

				#Min and Max
				if sts_switch[0] == True:
					for key in grid_data:
						if grid_data[key][z] > Fraction:
							grid_data[key][z] = Fraction

				if sts_switch[1] == True:
					if Fraction_Max[z] < Fraction:
						Fraction_Max[z] = Fraction

				#Total and Count for Mean
				if (sts_switch[2] == True) | (sts_switch[3] == True):
					TOT_Fraction[z] += Fraction
					Count[z] += 1 
				
				#Standard Deviation 
				if sts_switch[4] == True:
					TOT_Fraction_sq[z] += Fraction**2

				##1D Histogram 
				#if (sts_switch[5] == True) | (sts_switch[6] == True):
				#	hist_bnd1 = np.linspace(lobnd1,upbnd1,bin_num[0]+1)
				#	bin_interval1 = (upbnd1 - lobnd1)/bin_num[0]
				#	1D_hist_cnt = np.zeros(bin_num[0])
		#
				#	hist_idx1 = ((Fraction-lobnd1)/bin_interval1).astype(int)
				#	if hist_idx1 <= 1D_hist_cnt.shape[0]: 
				#		hist_idx1 = 1D_hist_cnt.shape[0]
				#	if hist_idx1 >= 0: 
				#		hist_idx1 = 0
				#	1D_hist_cnt[z, hist_idx1] += 1
#
				##2D Histogram 
				#if sts_switch[6] == True:
				#	hist_bnd2 = np.linspace(lobnd2,upbnd2,bin_num[1]+1)
				#	2D_hist_cnt = np.zeros((bin_num[0],bin_num[1]))
				#	bin_interval2 = (upbnd2 - lobnd2)/bin_num[1]
#
				#	hist_idx2 = ((Fraction-lobnd2)/bin_interval2).astype(int)
				#	if hist_idx2 <= hist_cnt2.shape[0]: 
				#		hist_idx2 = hist_cnt2.shape[0]
				#	if hist_idx2 >= 0: 
				#		hist_idx2 = 0
				#	2D_hist_cnt = [z, hist_idx1,hist_idx2] += 1


	return (Count,Fraction_Min,Fraction_Max,TOT_Fraction,TOT_Fraction_sq)


def addGridEntry(f,name,units,long_name,data):
	'''
	f:h5py.File()
	-------------------------------------
	Ex.
	self.addGridEntry(f,'CF','Fraction','Cloud_Fraction',total_cloud_fraction)
	'''
	PCentry=f.create_dataset(name,data=data)
	PCentry.dims[0].label='lat_bnd'
	PCentry.dims[1].label='lon_bnd'
	PCentry.attrs['units']=units
	PCentry.attrs["long_name"]=long_name	

if __name__ =='__main__':
# This is the main program for using concurrent to speed up the whole process
	
	#-------------STEP 0: Read the input from User --------
	# checking user input
	if (len(sys.argv) != 9) & (len(sys.argv) != 14):
		print("Wrong user input")
		print("usage: python aggre_stats_mpi.py <True/False> <True/False> <True/False> \
												<True/False> <True/False> <True/False> \
												<True/False> <Variable Name List> <Bin Size> \
												<JHist Variable Name> <lobnd1> <upbnd1> <lobnd2> <upbnd2>")
		sys.exit()
	else:
		# Define the statistics names for HDF5 output
		sts_name = ['Minimum','Maximum','Mean','Standard_Deviation', \
					'Pixel_Counts','Histogram_Counts','Jhisto_vs_']

		# Pass system arguments to the function
		sts_switch = np.array(sys.argv[1:8],dtype=np.int)
		sts_switch = np.array((sts_switch == 1))
		varlist  = sys.argv[8] 

		# Read the variable names from the variable name list
		text_file = open(varlist, "r")
		varnames  = text_file.read().split('\n')
			
		if (sts_switch[5] == True) | (sts_switch[6] == True):
			bin_num  = sys.argv[9] 
			hist_var = sys.argv[10]
			lobnd1   = sys.argv[11]
			upbnd1   = sys.argv[12]
			lobnd2   = sys.argv[13]
			upbnd2   = sys.argv[14]
			# Read the joint histogram names from the variable name list
			text_file = open(hist_var, "r")
			histnames = text_file.read().split('\n')
	

	#-------------STEP 1: Set up the specific directory --------
	MYD06_dir= '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD06_L2/'
	MYD06_prefix = 'MYD06_L2.A'
	MYD03_dir= '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD03/'
	MYD03_prefix = 'MYD03.A'
	fileformat = 'hdf'
	
	#-------------STEP 2: Set up spactial and temporal resolution & variable names----------
	NTA_lats = [-90,90]   #[  0,40] #[-90,90]   #[-30,30]    
	NTA_lons = [-180,180] #[-40,60] #[-180,180] #[-60,60]  
	
	gap_x, gap_y = 1,1 #0.625,0.5

	if ((NTA_lons[-1]-NTA_lons[0])%gap_x != 0) | ((NTA_lats[-1]-NTA_lats[0])%gap_y != 0): 
		print("Grid size should be dividable by the dimension of the selected region.")
		print("If you choose the region of latitude  from -40 to 40, then you gird size (gap_y) should be dividable by 80.")
		print("If you choose the region of longitude from  20 to 35, then you gird size (gap_x) should be dividable by 55.")
		print("Please try again!")
		sys.exit()

	map_lon = np.arange(NTA_lons[0],NTA_lons[1],gap_x)
	map_lat = np.arange(NTA_lats[0],NTA_lats[1],gap_y)
	Lon,Lat = np.meshgrid(map_lon,map_lat)
	
	grid_lon=np.int((NTA_lons[-1]-NTA_lons[0])/gap_x)
	grid_lat=np.int((NTA_lats[-1]-NTA_lats[0])/gap_y)

	#Count           = np.zeros(grid_lat*grid_lon)
	#Fraction_Min    = np.zeros(grid_lat*grid_lon) + np.inf
	#Fraction_Max    = np.zeros(grid_lat*grid_lon) - np.inf
	#TOT_Fraction    = np.zeros(grid_lat*grid_lon)
	#TOT_Fraction_sq = np.zeros(grid_lat*grid_lon)


	grid_data = {}
	sts_idx = np.array(np.where(sts_switch == True))[0]
	for key in varnames:
		for i in range(sts_idx.shape[0]):
			if i == 0:
				grid_data[key+'_'+sts_name[i]] = np.zeros(grid_lat*grid_lon) + np.inf
			elif i == 1:
				grid_data[key+'_'+sts_name[i]] = np.zeros(grid_lat*grid_lon) - np.inf
			elif i == 6:
				grid_data[histnames[0]+'_'+sts_name[i]+histnames[1]] = np.zeros(grid_lat*grid_lon)
			else:
				grid_data[key+'_'+sts_name[i]] = np.zeros(grid_lat*grid_lon)

	# Read the filename list for different time period
	fname1,fname2 = [],[]

	years  = np.array([2008])
	months = np.array([1])
	days = np.arange(1,2,dtype=np.int) 

	for yr,day in zip(years,days):
		yc ='%04i' % yr
		dc ='%03i' % day
		fname_tmp1 = read_filelist(MYD06_dir,MYD06_prefix,yc,dc,fileformat)
		fname_tmp2 = read_filelist(MYD03_dir,MYD03_prefix,yc,dc,fileformat)
		fname1 = np.append(fname1,fname_tmp1)
		fname2 = np.append(fname2,fname_tmp2)

	filenum = np.arange(len(fname1))

	# Start counting operation time
	start_time = timeit.default_timer() 

	results = np.asarray(run_modis_aggre(fname1,fname2,NTA_lats,NTA_lons,grid_lon,gap_x,gap_y,filenum, \
										 grid_data,sts_switch,varnames,bin_num,lobnd1,upbnd1,lobnd2,upbnd2))
		
	Count           += results[0,:]
	Fraction_Min     = results[1,:]
	Fraction_Max     = results[2,:]
	TOT_Fraction    += results[3,:]
	TOT_Fraction_sq += results[4,:]

	for i in range(1,size):
		recv_req = comm.Irecv(results,source=i, tag=0)
		recv_req.wait()
		
		Count = Count + results[1,:]
		Fraction_Min = np.dstack((Fraction_Min,results[2,:]))
		Fraction_Max = np.dstack((Fraction_Max,results[3,:]))
		TOT_Fraction = TOT_Fraction + results[0,:]
		TOT_Fraction_sq = TOT_Fraction_sq + results[4,:]

	# Compute the mean cloud fraction & Statistics (Include Min & Max & Standard deviation)
	Mean_Fraction = (TOT_Fraction / Count)
	Std_Fraction  = (TOT_Fraction_sq / Count) - Mean_Fraction**2

	Count         =         Count.reshape([grid_lat,grid_lon])
	Mean_Fraction = Mean_Fraction.reshape([grid_lat,grid_lon])
	Std_Fraction  =  Std_Fraction.reshape([grid_lat,grid_lon])

	Fraction_Min = np.min(Fraction_Min,axis=2).reshape([grid_lat,grid_lon])
	Fraction_Max = np.max(Fraction_Max,axis=2).reshape([grid_lat,grid_lon])

	end_time = timeit.default_timer()

	print('Mean_Fraction:')
	print( Mean_Fraction  )

	print ("Operation Time in {:7.2f} seconds".format(end_time - start_time))
	
	# Create HDF5 file to store the result 
	l3name='MOD08_M3'+'A{:04d}{:02d}'.format(years[0],months[0])
	ff=h5py.File(l3name+'.hdf5','w')

	PC=ff.create_dataset('lat_bnd',data=map_lat)
	PC.attrs['units']='degrees'
	PC.attrs['long_name']='Latitude_boundaries'    

	PC=ff.create_dataset('lon_bnd',data=map_lon)
	PC.attrs['units']='degrees'
	PC.attrs['long_name']='Longitude_boundaries'    
	
	addGridEntry(ff,'Cloud_Fraction_Mean'              ,'none','Cloud Fraction from Cloud Mask (cloudy & prob cloudy)',Mean_Fraction)
	addGridEntry(ff,'Cloud_Fraction_Standard_Deviation','none','Cloud Fraction from Cloud Mask (cloudy & prob cloudy)',Std_Fraction )
	addGridEntry(ff,'Cloud_Fraction_Minimum'           ,'none','Cloud Fraction from Cloud Mask (cloudy & prob cloudy)',Fraction_Min )
	addGridEntry(ff,'Cloud_Fraction_Maximum'           ,'none','Cloud Fraction from Cloud Mask (cloudy & prob cloudy)',Fraction_Max )
	addGridEntry(ff,'Cloud_Fraction_Pixel_Counts'      ,'none','Cloud Fraction from Cloud Mask (cloudy & prob cloudy)',Count) 
	
	ff.close()

	print(l3name+'.hdf5 Saved!')
