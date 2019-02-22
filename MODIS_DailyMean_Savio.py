#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 10:16:55 2019
@author: saviokay
"""
import time
import math
import os,datetime,sys,fnmatch
from netCDF4 import Dataset
import numpy as np
import xarray as xr
import pandas as pd
import dask.array as da

def read_MODIS_level2_data(MOD06_file,MOD03_file):
    print(MOD06_file)
    print(MOD03_file)
    print('reading the cloud mask from MOD06_L2 product')
    myd06 = Dataset(MOD06_file, "r")
    CM = myd06.variables["Cloud_Mask_1km"][:,:,:] # Reading Specific Variable 'Cloud_Mask_1km'.
    CM   = (np.array(CM[:,:,0],dtype='byte') & 0b00000110) >>1
    CM = np.array(CM).byteswap().newbyteorder()
    #CM = np.ravel(CM) # Changing The Shape Of The Variable 'Longitude'.
    print('level-2 cloud mask array shape',CM.shape)

    myd03 = Dataset(MOD03_file, "r")
    print('reading the lat-lon from MOD03 product')
    latitude = myd03.variables["Latitude"][:,:] # Reading Specific Variable 'Latitude'.
    #latitude = np.ravel(latitude) # Changing The Shape Of The Variable 'Latitude'.
    latitude = np.array(latitude).byteswap().newbyteorder() # Addressing Byteswap For Big Endian Error.
    longitude = myd03.variables["Longitude"][:,:] # Reading Specific Variable 'Longitude'.
    #longitude = np.ravel(longitude) # Changing The Shape Of The Variable 'Longitude'.
    longitude = np.array(longitude).byteswap().newbyteorder() # Addressing Byteswap For Big Endian Error.
    print('level-2 lat-lon array shape',latitude.shape)

    return latitude,longitude,CM

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


if __name__ == '__main__':
    import itertools
    MOD03_path = 'HDFFiles/'
    MOD06_path = 'HDFFiles/'
    satellite = 'Aqua'

    yr = [2008]
    mn = [1] #np.arange(1,13)  #[1]
    dy = [1] #np.arange(1,32) # [1]
    #np.arange(1,31)
    # latitude and longtitude boundaries of level-3 grid
    lat_bnd = np.arange(-90,91,1)
    lon_bnd = np.arange(-180,180,1)
    nlat = 180
    nlon = 360

    TOT_pix      = np.zeros(nlat*nlon)
    CLD_pix      = np.zeros(nlat*nlon)

    MOD03_fp = 'MYD03.A*.hdf'
    MOD06_fp = 'MYD06_L2.A*.hdf'
    MOD03_fn, MOD06_fn =[],[]
    #MOD03_fn2, MOD06_fn2 =[],[]
    MOD03_fn2, MOD06_fn2 =[],[]
    for MOD06_flist in  os.listdir(MOD06_path):
        if fnmatch.fnmatch(MOD06_flist, MOD06_fp):
            MOD06_fn = MOD06_flist
            MOD06_fn2.append(MOD06_flist)
            print(MOD06_fn)
    for MOD03_flist in  os.listdir(MOD03_path):
        if fnmatch.fnmatch(MOD03_flist, MOD03_fp):
            MOD03_fn = MOD03_flist
            MOD03_fn2.append(MOD03_flist)
            print(MOD03_fn)
    if MOD03_fn and MOD06_fn:
        # if both MOD06 and MOD03 products are in the directory
        print('reading level 2 geolocation and cloud data')
        #print(MOD06_fn)
        #print(MOD03_fn)
        Lat,Lon,CM = read_MODIS_level2_data(MOD06_path+MOD06_fn,MOD03_path+MOD03_fn)



    myd06_name = MOD06_path
    #cm = np.array(2)
    cm = np.zeros((2030,1354), dtype=np.float32)
    #cm = np.empty([2, 2])
    #cm1 = np.array(2)
    cm1 = np.zeros((2030,1354), dtype=np.float32)
    #cm1 = np.empty([2, 2])
    #cmr = np.array(2)
    cmr = np.zeros((2030,1354), dtype=np.float32)
    #cmr = np.empty([2, 2])
    for MOD06_file in MOD06_fn2:
        MOD06_file2 = myd06_name + MOD06_file
        myd06 = Dataset(MOD06_file2, "r")
        CM = myd06.variables["Cloud_Mask_1km"][:,:,:]# Reading Specific Variable 'Cloud_Mask_1km'.
        CMr = myd06.variables["Cloud_Mask_1km"][:,:,0]
        CM1   = (np.array(CM[:,:,0],dtype='byte') & 0b00000001)
        CM   = (np.array(CM[:,:,0],dtype='byte') & 0b00000110) >>1
        CM = np.array(CM).byteswap().newbyteorder()
        cm = np.dstack((cm,CM))
        cm1 = np.dstack((cm1,CM1))
        cmr = np.dstack((cmr,CMr))
        #CM = np.ravel(CM) # Changing The Shape Of The Variable 'Longitude'.
        #ar.append(CM)
        #cm = np.append(cm, CM)
        #cm = np.concatenate((cm,CM), axis=0)
        #cm = np.column_stack([cm,CM])
        #cm1 = np.append(cm1, CM1)
        #cmr = np.append(cmr, CMr)
        #ar.insert(CM)
        print('level-2 cloud mask array shape',CM.shape)


    myd03_name = MOD03_path
    #lat = np.array(2)
    #lon = np.array(2)
    lat = np.zeros((2030,1354), dtype=np.float32)
    lon = np.zeros((2030,1354), dtype=np.float32)
    for MOD03_file in MOD03_fn2:
        MOD03_file2 = myd03_name + MOD03_file
        myd03 = Dataset(MOD03_file2, "r")
        latitude = myd03.variables["Latitude"][:,:] # Reading Specific Variable 'Latitude'.
        lat = np.dstack((lat,latitude))
        #latitude = np.ravel(latitude) # Changing The Shape Of The Variable 'Latitude'.
        #latitude = np.array(latitude).byteswap().newbyteorder() # Addressing Byteswap For Big Endian Error.
        #lat = np.append(lat, latitude)
        print('Latitude Shape Is: ',lat.shape)

        longitude = myd03.variables["Longitude"][:,:] # Reading Specific Variable 'Longitude'.
        #longitude = np.ravel(longitude) # Changing The Shape Of The Variable 'Longitude'.
        #longitude = np.array(longitude).byteswap().newbyteorder() # Addressing Byteswap For Big Endian Error.
        lon = np.dstack((lon,longitude))
        print('Latitude Shape Is: ',lon.shape)


    mod03_num_chunk = len(MOD03_fn2)
    lat = da.from_array(lat, chunks =(2030,1354,1))
    lon = da.from_array(lon, chunks =(2030,1354,1))
    CMB123 = lat.astype(np.int8)
    CMB124 = lon.astype(np.int8)
    #CMB123 = xr.DataArray(CMB123)
    #CMB124 = xr.DataArray(CMB124)
    mod06_num_chunk = len(MOD06_fn2)
    cm = da.from_array(cm, chunks =(2030,1354,1))
    cm1 = da.from_array(cm1, chunks =(2030,1354,1))
    cmr = da.from_array(cmr, chunks =(2030,1354,1))
    CMB123 = da.from_array(CMB123, chunks =(2030,1354,1))
    CMB124 = da.from_array(CMB124, chunks =(2030,1354,1))
    CMrd = cmr.compute()


    #--- Creating Xarray ---#
    dsa3 = xr.Dataset()
    dsa3.coords['Latitude'] = (('x','y','z'), lat)
    dsa3.coords['Longitude'] = (('x','y','z'), lon)
    dsa3.coords['Bit0'] = (('x','y','z'), cm1)
    dsa3.coords['Bit12'] = (('x','y','z'), cm)
    dsa3.coords['LatInt'] = (('x','y','z'), CMB123)
    dsa3.coords['LonInt'] = (('x','y','z'), CMB124)
    dsa3['CloudFraction'] = (('x','y','z'), CMrd)

    Nobs = np.sum(cm1.compute()>=0)
    print("Number Of Total Pixel: ", np.sum(cm1.compute()>=0), " Pixels")
    print("Number Of Pixel With Cloud Mask Determined: ", np.sum(cm1.compute()==1), " Pixels")

    dsa4 = list(dsa3.groupby('LatInt', 'LonInt'))

    list1=[]
    for key, value in dsa4:
        #print(key)
        list1.append(key)

    from dask.distributed import Client
    client = Client("tcp://10.49.83.13:41751")
    client

    i = 0
    start_time = time.time()
    while i < len(list1):
        #print(list1[i])
        a = list1[i]
        #print(dsa4[i])
        dsa4 =  list(dsa3.groupby('LatInt', 'LonInt'))
        dsa41 = dsa4[i][1]
        dsa51 = list(dsa41.groupby('LonInt'))
        dsa61 = dsa51[0][1]
        dsa71 = list(dsa61.groupby('Bit0', 'Bit12'))
        dsa81 = dsa71[0][1]
        dsa91 = list(dsa81.groupby('Bit12'))
        dsa81 = dsa71[0][1]
        dsa91 = dsa81.groupby('Bit12').count()
        print(dsa91)
        i += 1
    end_time = time.time()
    print("Total Time Taken This Loop: ", end_time - start_time)
    hours, rem = divmod(end_time-start_time, 3600)
    minutes, seconds = divmod(rem, 60)
    print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))
