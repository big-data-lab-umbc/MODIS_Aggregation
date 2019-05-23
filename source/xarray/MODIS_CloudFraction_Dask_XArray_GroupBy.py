#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 04 15:26:55 2019
@author: jianwuwang, saviokay
"""
import numpy as np
import xarray as xr
import pandas as pd
from netCDF4 import Dataset
import os,datetime,sys,fnmatch
import time
import math
import dask.array as da
import h5py
from dask_jobqueue import SLURMCluster
from dask.distributed import Client
import argparse


#Function: Read Function For MODO3 & MODO6 Files.
def read_MODIS_level2_data(MOD06_file,MOD03_file):
    #print('Reading The Cloud Mask From MOD06_L2 Product:')
    myd06 = Dataset(MOD06_file, "r")
    CM = myd06.variables["Cloud_Mask_1km"][:,:,:] # Reading Specific Variable 'Cloud_Mask_1km'.
    CM   = (np.array(CM[:,:,0],dtype='byte') & 0b00000110) >>1
    CM = np.array(CM).byteswap().newbyteorder()
   # print('The Level-2 Cloud Mask Array Shape',CM.shape)
    print(' ')

    myd03 = Dataset(MOD03_file, "r")
    #print('Reading The Latitude-Longitude From MOD03 Product:')
    latitude = myd03.variables["Latitude"][:,:] # Reading Specific Variable 'Latitude'.
    latitude = np.array(latitude).byteswap().newbyteorder() # Addressing Byteswap For Big Endian Error.
    longitude = myd03.variables["Longitude"][:,:] # Reading Specific Variable 'Longitude'.
    longitude = np.array(longitude).byteswap().newbyteorder() # Addressing Byteswap For Big Endian Error.
    #print('The Level-2 Latitude-Longitude Array Shape',latitude.shape)
    print(' ')

    return latitude,longitude,CM

#Function: Filtering Function For Cloud Fraction.
def filter_cf(x, axis=1):
    count0 = 0
    for i in x:
        if i <= 1:
            count0 +=1
    return count0/len(x)

#Function: Save Total Cloud Fraction And GeoLocation Varibales To HDF5 Output File.
def save_hdf(out_name,total_cloud_fraction,lat_bnd,lon_bnd):
    f=h5py.File(out_name,'w')
    PCentry=f.create_dataset('CF',data=total_cloud_fraction)
    PCentry.dims[0].label='lat_bnd'
    PCentry.dims[1].label='lon_bnd'

    PC=f.create_dataset('lat_bnd',data=lat_bnd)
    PC.attrs['units']='degrees'
    PC.attrs['long_name']='Latitude_boundaries'

    PC=f.create_dataset('lon_bnd',data=lon_bnd)
    PC.attrs['units']='degrees'
    PC.attrs['long_name']='Longitude_boundaries'
    f.close()
    print(out_name+' Saved!!')


# Beginning of the main code

# Read input from commmandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('mod3_path', nargs='?', const=1, default="input-data/MYD03/")
parser.add_argument('mod6_path', nargs='?', const=1, default="input-data/MYD06/")
parser.add_argument('output_path', nargs='?', const=1, default="output-data/output.hdf5")
args = parser.parse_args()
MOD03_path =args.mod3_path #'input-data/MYD3'
MOD06_path =args.mod6_path #'input-data/MYD6'
outfile_name = args.output_path #'output-data/'


satellite = 'Aqua'
yr = [2008]
mn = [1]
dy = [1]
# Latitude & Longtitude Boundaries Of Level-3 Grid
lat_bnd = np.arange(-90,91,1)
lon_bnd = np.arange(-180,180,1)
nlat = 180
nlon = 360

MOD03_filepath = 'MYD03.A*.hdf'
MOD06_filepath = 'MYD06_L2.A*.hdf'
MOD03_filename, MOD06_filename =[],[]
MOD03_filename2, MOD06_filename2 =[],[]

for MOD06_filelist in  os.listdir(MOD06_path):
    if fnmatch.fnmatch(MOD06_filelist, MOD06_filepath):
        MOD06_filename = MOD06_filelist
        MOD06_filename2.append(MOD06_filelist)

for MOD03_filelist in  os.listdir(MOD03_path):
    if fnmatch.fnmatch(MOD03_filelist, MOD03_filepath):
        MOD03_filename = MOD03_filelist
        MOD03_filename2.append(MOD03_filelist)

#if MOD03_filename and MOD06_filename:
#    print('Reading Level 2 GeoLocation & Cloud Data')
#    Lat,Lon,CM = read_MODIS_level2_data(MOD06_path+MOD06_filename,MOD03_path+MOD03_filename)

print('The Number Of Files In The MODO3 List: ')
print(len(MOD03_filename2))
print(' ')
print('The Number Of Files In The MODO6_L2 List: ')
print(len(MOD06_filename2))
print(' ')

cm = np.zeros((2030,1354), dtype=np.float32)

for MOD06_file in MOD06_filename2:
    MOD06_file2 = MOD06_path + MOD06_file
    myd06 = Dataset(MOD06_file2, "r")
    # Reading Specific Variable 'Cloud_Mask_1km'.
    CM = myd06.variables["Cloud_Mask_1km"][:,:,0]
    CM = (np.array(CM,dtype='byte') & 0b00000110) >>1
    CM = np.array(CM).byteswap().newbyteorder()
    cm = np.dstack((cm,CM))

print('The Cloud Mask Array Shape Is: ',cm.shape)

lat = np.zeros((2030,1354), dtype=np.float32)
lon = np.zeros((2030,1354), dtype=np.float32)

for MOD03_file in MOD03_filename2:
    MOD03_file2 = MOD03_path + MOD03_file
    myd03 = Dataset(MOD03_file2, "r")
    latitude = myd03.variables["Latitude"][:,:] # Reading Specific Variable 'Latitude'.
    lat = np.dstack((lat,latitude))

    longitude = myd03.variables["Longitude"][:,:] # Reading Specific Variable 'Longitude'.
    lon = np.dstack((lon,longitude))



print('The Latitude Array Shape Is: ',lat.shape)
print('The Longitude Array Shape Is: ',lon.shape)

# Setting Integer Values Of Latitude & Longitude For GroupBy Operation.
latint = lat.astype(np.int8)
lonint = lon.astype(np.int8)

# Creating Cluster With The Function 'SLURMCluster' With Cluster's Configurations As Arguments.
#cluster = SLURMCluster(cores=16, memory='32GB', project='pi_jianwu', walltime='02:00:00', queue='batch', job_extra=['--qos=medium+'])
#cluster.scale(5) #Scaling To With Scale Function.

#client = Client(cluster)
#print(client)
#print(cluster)

# Creating Dask Arrays With Chunks With The Array Variables.
lat = da.from_array(lat, chunks =(2030,1354,1))
lon = da.from_array(lon, chunks =(2030,1354,1))
latint = da.from_array(latint, chunks =(2030,1354,1))
lonint = da.from_array(lonint, chunks =(2030,1354,1))
cm = da.from_array(cm, chunks =(2030,1354,1))

# Creating XArrays Data Structure With CloudMask As The Data Variable.
dsa = xr.Dataset()
dsa.coords['Latitude'] = (('x','y','z'), lat)
dsa.coords['Longitude'] = (('x','y','z'), lon)
dsa.coords['LatInt'] = (('x','y','z'), latint)
dsa.coords['LonInt'] = (('x','y','z'), lonint)
dsa['CloudMask'] = (('x','y','z'), cm)

# Performing GroupBy Operation To Obtain The Cloud Fractions.
start_time = time.time()
finallist = [] #Final List Of Latitude-Longitude-CM
dsag = list(dsa.groupby('LonInt'))
for listOfOneLong in dsag:
    oneGrid = listOfOneLong[1].groupby('LatInt').reduce(filter_cf)
    for index in range(0, oneGrid.CloudMask.size-1):
        finalvalues = [oneGrid.LatInt.data[index], listOfOneLong[0], oneGrid.CloudMask.data[index]]
        finallist.append(finalvalues)

end_time = time.time()
print("Total Time Taken This Loop: ", end_time - start_time)
hours, rem = divmod(end_time-start_time, 3600)
minutes, seconds = divmod(rem, 60)
print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))

cf_array = np.zeros((180,360))

# Creating Final List Of Latitude-Longitude-CM.
cf_array[:]=np.nan
print('The Length Of The Final List Of Latitude-Longitude-CM: ')
print(len(finallist))
print(' ')
print('-------------------------------------')
for i in range(len(finallist)):
    cf_array[(90-finallist[i][0]),(finallist[i][1]+180)] = finallist[i][2]

print(cf_array)

# Read Total Cloud Fraction And GeoLocation Varibales To Create/Save HDF5 Output File.
total_cloud_fraction = cf_array
save_hdf(outfile_name,total_cloud_fraction,lat_bnd,lon_bnd)

# Final Maintenance.
del dsag
del oneGrid
#cluster.close()


