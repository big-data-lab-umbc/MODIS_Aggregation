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

def read_filelist(loc_dir,prefix,yr,day,fileformat):
	# Read the filelist in the specific directory
	str = os.popen("ls "+ loc_dir + prefix + yr + day + "*."+fileformat).read()
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

def run_modis_aggre(fname1,fname2,NTA_lats,NTA_lons,grid_lon,gap_x,gap_y,hdfs):
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
				if Fraction_Min[z] > Fraction:
					Fraction_Min[z] = Fraction
				if Fraction_Max[z] < Fraction:
					Fraction_Max[z] = Fraction

				#Total and Count for Mean
				TOT_Fraction[z] += Fraction
				Count[z] += 1 
				
				#Standard Deviation 
				TOT_Fraction_sq[z] += Fraction**2

	return (Count,Fraction_Min,Fraction_Max,TOT_Fraction,TOT_Fraction_sq)

#Mean and std. computations. minmax() is called inside this function
def MeanStd(data,z,latlon_index,M):
	#Both mean and stdd
	#print(key)
	val=data[np.where(latlon_index == z)]
	M.XXX_pix[key][z]=M.XXX_pix[key][z]+np.sum(val)
	M.XXX_pixSq[key][z]=M.XXX_pixSq[key][z]+np.sum(val**2)   
	minmax(val,z,M)

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
	
	#-------------STEP 1: Set up the specific directory --------
	MYD06_dir= '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD06_L2/'
	MYD06_prefix = 'MYD06_L2.A'
	MYD03_dir= '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD03/'
	MYD03_prefix = 'MYD03.A'
	fileformat = 'hdf'
	
	#-------------STEP 2: Set up spactial and temporal resolution----------
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

	#print(grid_lon,grid_lat,grid_lat*grid_lon)
	
	Count           = np.zeros(grid_lat*grid_lon)
	Fraction_Min    = np.zeros(grid_lat*grid_lon) + np.inf
	Fraction_Max    = np.zeros(grid_lat*grid_lon) - np.inf
	TOT_Fraction    = np.zeros(grid_lat*grid_lon)
	TOT_Fraction_sq = np.zeros(grid_lat*grid_lon)
	
	fname1,fname2 = [],[]

	# Read all files in a month (in this case: January)
	# Read the filename list for different time period
	#days = np.arange(1,2,dtype=np.int)
	years  = np.array([2008])
	months = np.array([1])
	days = np.arange(1,2,dtype=np.int) 

	for yr,day in days,years:
		yc ='%04i' % iy
		dc ='%03i' % day
		fname_tmp1 = read_filelist(MYD06_dir,MYD06_prefix,yr,dc,fileformat)
		fname_tmp2 = read_filelist(MYD03_dir,MYD03_prefix,yr,dc,fileformat)
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

	results = np.asarray(run_modis_aggre(fname1,fname2,NTA_lats,NTA_lons,grid_lon,gap_x,gap_y,hdfs))
		
	if rank == 0:
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
		
		# Create file to store the result 
		#np.savetxt("cloud_fraction_mean.dat", Mean_Fraction, fmt="%10.4f")
		#np.savetxt("cloud_fraction_min.dat" , Fraction_Min , fmt="%10.4f")
		#np.savetxt("cloud_fraction_max.dat" , Fraction_Max , fmt="%10.4f")
		#np.savetxt("cloud_fraction_std.dat" , Std_Fraction , fmt="%10.4f")
		#np.savetxt("cloud_fraction_pix_count.dat",   Count , fmt="%10d")
		#np.savetxt("test_geolocation_lat.dat" , Lat, fmt="%10.4f")
		#np.savetxt("test_geolocation_lon.dat" , Lon, fmt="%10.4f")

		# Create HDF5 file to store the result 
		l3name='MOD08_M3'+'A{:04d}{:02d}'.format(years[0],months[0])
		ff=h5py.File(l3name+'.hdf5','w')
		addGridEntry(ff,'Cloud_Fraction_Mean'              ,'none','Cloud Fraction from Cloud Mask (cloudy & prob cloudy)',Mean_Fraction)
		addGridEntry(ff,'Cloud_Fraction_Standard_Deviation','none','Cloud Fraction from Cloud Mask (cloudy & prob cloudy)',Std_Fraction )
		addGridEntry(ff,'Cloud_Fraction_Minimum'           ,'none','Cloud Fraction from Cloud Mask (cloudy & prob cloudy)',Fraction_Min )
		addGridEntry(ff,'Cloud_Fraction_Maximum'           ,'none','Cloud Fraction from Cloud Mask (cloudy & prob cloudy)',Fraction_Max )
		addGridEntry(ff,'Cloud_Fraction_Pixel_Counts'      ,'none','Cloud Fraction from Cloud Mask (cloudy & prob cloudy)',Count) 
		
		PC=ff.create_dataset('lat_bnd',data=Agg.map_lat)
        PC.attrs['units']='degrees'
        PC.attrs['long_name']='Latitude_boundaries'    

        PC=ff.create_dataset('lon_bnd',data=Agg.map_lon)
        PC.attrs['units']='degrees'
        PC.attrs['long_name']='Longitude_boundaries'    
        ff.close()

        print(l3name+'.hdf5 Saved!')

def save_level3_hdf5(self,Agg):
        '''
        To save aggregated data products.
        Agg: MODIS_L2toL3 object
        '''
        self.MODIS_L2toL3=Agg
        self.fname=Agg.l3name
        ff=h5py.File(self.fname+'.hdf5','w')
        self.addGridEntry(ff,'CF','Fraction','Cloud_Fraction',Agg.M.total_cloud_fraction)
        self.addGridEntry(ff,'PC','Count','Pixel_Count',Agg.M.pixel_count)
        for key in Agg.variables:
            for st in Agg.M.stt:
                self.addGridEntry(ff, key+'_'+st, Agg.variables[key][1], Agg.variables[key][0]+'_'+self.get_long_name(st), \
                                  Agg.M.stt[st][key])
        PC=ff.create_dataset('lat_bnd',data=Agg.lat_bnd)
        PC.attrs['units']='degrees'
        PC.attrs['long_name']='Latitude_boundaries'    

        PC=ff.create_dataset('lon_bnd',data=Agg.lon_bnd)
        PC.attrs['units']='degrees'
        PC.attrs['long_name']='Longitude_boundaries'    
        ff.close()
        print(self.fname+'.hdf5 Saved!')



	else:
		print("Process {} finished".format(rank))
		send_req = comm.Isend(results, dest=0, tag=0)
		send_req.wait()