#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  22 17:16:55 2017
@author: saviokay
"""
import numpy as np
import xarray as xr
import pandas as pd
from netCDF4 import Dataset
import os,datetime,sys,fnmatch
import time
import math
import dask.array as da

from dask_jobqueue import SLURMCluster
cluster = SLURMCluster(cores=16, memory='32GB', project='pi_jianwu', walltime='02:00:00', queue='batch', job_extra=['--qos=medium+'])
cluster.scale(5)
from dask.distributed import Client
client = Client(cluster)
print(client)
print(cluster)


# Read Function For MODO3 & MODO6 Files.
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

# Misc Function For Processing Cloud Fraction.
def value_locate(refx, x):
    refx = np.array(refx)
    x = np.array(x)
    loc = np.zeros(len(x), dtype='int')

    for i in range(len(x)):
        ix = x[i]
        ind = ((refx - ix) <= 0).nonzero()[0]
        if len(ind) == 0:
            loc[i] = -1
        else: loc[i] = ind[-1]

    return loc

def division(n, d):

    div = np.zeros(len(d))
    for i in range(len(d)):
        if d[i] >0:
          div[i]=n[i]/d[i]
        else: div[i]=None 

    return div

def countzero(x, axis=1):
    #print(x)
    count0 = 0
    count1 = 0
    for i in x:
        if i <= 1:
            count0 +=1
    #print(count0/len(x))
    return count0/len(x)

# Setting File Location As Environment Variables
#MOD03_path = os.environ.get('MOD03_PATH')
#MOD06_path = os.environ.get('MOD06_PATH')

MOD03_path =sys.argv[1] #'HDFFiles/'
MOD06_path =sys.argv[2] #'HDFFiles/'


satellite = 'Aqua'

yr = [2008]
mn = [1] #np.arange(1,13)  #[1]
dy = [1] #np.arange(1,32) # [1] #np.arange(1,31)
lat_bnd = np.arange(-90,91,1)# latitude and longtitude boundaries of level-3 grid
lon_bnd = np.arange(-180,180,1)
nlat = 180
nlon = 360

MOD03_fp = 'MYD03.A*.hdf'
MOD06_fp = 'MYD06_L2.A*.hdf'
MOD03_fn, MOD06_fn =[],[]
#MOD03_fn2, MOD06_fn2 =[],[]
MOD03_fn2, MOD06_fn2 =[],[]
for MOD06_flist in  os.listdir(MOD06_path):
    if fnmatch.fnmatch(MOD06_flist, MOD06_fp):
        MOD06_fn = MOD06_flist
        MOD06_fn2.append(MOD06_flist)
   
for MOD03_flist in  os.listdir(MOD03_path):
    if fnmatch.fnmatch(MOD03_flist, MOD03_fp):
        MOD03_fn = MOD03_flist
        MOD03_fn2.append(MOD03_flist)

if MOD03_fn and MOD06_fn:
    print('Reading Level 2 GeoLocation & Cloud Data')   
    Lat,Lon,CM = read_MODIS_level2_data(MOD06_path+MOD06_fn,MOD03_path+MOD03_fn)

print('The Number Of Files In The MODO3 List: ')
print(len(MOD03_fn2))
print(' ')
print('The Number Of Files In The MODO6_L2 List: ')
print(len(MOD06_fn2))
print(' ')

myd06_name = sys.argv[1]#'HDFFiles/'
cm = np.zeros((2030,1354), dtype=np.float32)
#cm = np.empty((2030,1354), dtype=np.float32)
bit0 = np.zeros((2030,1354), dtype=np.float32)
bit12 = np.zeros((2030,1354), dtype=np.float32)

for MOD06_file in MOD06_fn2:
    MOD06_file2 = myd06_name + MOD06_file
    myd06 = Dataset(MOD06_file2, "r")
    CM = myd06.variables["Cloud_Mask_1km"][:,:,:]# Reading Specific Variable 'Cloud_Mask_1km'.
    CM = myd06.variables["Cloud_Mask_1km"][:,:,0]
    bit0r   = (np.array(CM,dtype='byte') & 0b00000001)
    bit12r   = (np.array(CM,dtype='byte') & 0b00000110) >>1
    CM = np.array(CM).byteswap().newbyteorder()
    cm = np.dstack((cm,CM))
    bit0 = np.dstack((bit0,bit0r))
    bit12 = np.dstack((bit12,bit12r))
        
print('The Cloud Mask Array Shape Is: ',cm.shape)

#myd03_name = os.environ.get('MOD03_PATH')
myd03_name = sys.argv[2] #'HDFFiles/'


lat = np.zeros((2030,1354), dtype=np.float32)
lon = np.zeros((2030,1354), dtype=np.float32)

for MOD03_file in MOD03_fn2:
    MOD03_file2 = myd03_name + MOD03_file
    myd03 = Dataset(MOD03_file2, "r")
    latitude = myd03.variables["Latitude"][:,:] # Reading Specific Variable 'Latitude'.
    lat = np.dstack((lat,latitude))

    longitude = myd03.variables["Longitude"][:,:] # Reading Specific Variable 'Longitude'.
    lon = np.dstack((lon,longitude))


    
print('The Latitude Array Shape Is: ',lat.shape)
print('The Longitude Array Shape Is: ',lon.shape)

latint = lat.astype(np.int8)
lonint = lon.astype(np.int8)

lat = da.from_array(lat, chunks =(2030,1354,1))
lon = da.from_array(lon, chunks =(2030,1354,1))
latint = da.from_array(latint, chunks =(2030,1354,1))
lonint = da.from_array(lonint, chunks =(2030,1354,1))
cm = da.from_array(cm, chunks =(2030,1354,1))
bit0 = da.from_array(bit0, chunks =(2030,1354,1))
bit12 = da.from_array(bit12, chunks =(2030,1354,1))

dsa = xr.Dataset()
dsa.coords['Latitude'] = (('x','y','z'), lat)
dsa.coords['Longitude'] = (('x','y','z'), lon)
dsa.coords['LatInt'] = (('x','y','z'), latint)
dsa.coords['LonInt'] = (('x','y','z'), lonint)
dsa['CloudMask'] = (('x','y','z'), cm)

'''
from dask_jobqueue import SLURMCluster
cluster = SLURMCluster(cores=16, memory='32GB', project='pi_jianwu', walltime='02:00:00', queue='batch', job_extra=['--qos=medium+'])
cluster.scale(5)
from dask.distributed import Client
client = Client(cluster)
client
'''
start_time = time.time()

listOfCM = []
finallist = []


dsag = list(dsa.groupby('LonInt'))
for listOfOneLong in dsag:
    oneGrid = listOfOneLong[1].groupby('LatInt').reduce(countzero)
    for index in range(0, oneGrid.CloudMask.size-1):
        finalvalues = [oneGrid.LatInt.data[index], listOfOneLong[0], oneGrid.CloudMask.data[index]]
        finallist.append(finalvalues)
       
end_time = time.time()
print("Total Time Taken This Loop: ", end_time - start_time)
hours, rem = divmod(end_time-start_time, 3600)
minutes, seconds = divmod(rem, 60)
print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))
lat_array = np.arange(-90,91,1)
lon_array = np.arange(-180,180,1)
cf_array = np.zeros((180,360))

cf_array[:]=np.nan
print('The Length Of The Final List Of Latitude-Longitude-CM: ')
print(len(finallist))
print(' ')
print('-------------------------------------')
for i in range(len(finallist)):
    cf_array[(90-finallist[i][0]),(finallist[i][1]+180)] = finallist[i][2]

print(cf_array)
#print(np.count_nonzero(~np.isnan(cf_array)))

import h5py
total_cloud_fraction = cf_array
out_name = 'output_final3.hdf5'
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

save_hdf(out_name,total_cloud_fraction,lat_bnd,lon_bnd)
cluster.close()
