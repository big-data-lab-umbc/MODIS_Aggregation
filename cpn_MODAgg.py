#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Sun Apr 21 13:38:21 2019
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
MODIS data aggregation - Multiple variables
"""
from netCDF4 import Dataset
import numpy as np
import h5py

def readEntry(key,ncf):
    ncf.variables[key][:]
    rdval=ncf.variables[key][:]
    scale=ncf.variables[key].getncattr('scale_factor')
    offst=ncf.variables[key].getncattr('add_offset')
    return (rdval+offst)*scale

def read_MODIS_level2_dataV2(MYD06_file,MYD03_file):
    #Reading variables from MYD06
    myd06 = Dataset(MYD06_file, "r")
    CM1km  = readEntry('Cloud_Mask_1km',myd06)             #Cloud mask
    CM     = (np.array(CM1km[:,:,0],dtype='byte') & 0b00000110) >>1
    CTP = readEntry('cloud_top_pressure_1km',myd06)     #Cloud Top Pressure (hPa)
    CTT = readEntry('cloud_top_temperature_1km',myd06)  #Cloud Top Temperature (K)
    CTH = readEntry('cloud_top_height_1km',myd06)       #Cloud Top Height (m)
    
    myd03 = Dataset(MYD03_file, "r")
    latitude = myd03.variables["Latitude"][:,:] # Reading Specific Variable 'Latitude'.
    latitude = np.array(latitude).byteswap().newbyteorder() # Addressing Byteswap For Big Endian Error.
    longitude = myd03.variables["Longitude"][:,:] # Reading Specific Variable 'Longitude'.
    longitude = np.array(longitude).byteswap().newbyteorder() # Addressing Byteswap For Big Endian Error.
    data={'CM':CM,'CTP':CTP,'CTT':CTT,'CTH':CTH}
    return latitude,longitude,data
def save_hdf(out_name,total_cloud_fraction,mean,lat_bnd,lon_bnd):
    f=h5py.File(out_name,'w')
    PCentry=f.create_dataset('CF',data=total_cloud_fraction)
    PCentry.dims[0].label='lat_bnd'
    PCentry.dims[1].label='lon_bnd'
    
    PCentry=f.create_dataset('CTP',data=mean['CTP'])
    PCentry.dims[0].label='lat_bnd'
    PCentry.dims[1].label='lon_bnd'
    
    PCentry=f.create_dataset('CTT',data=mean['CTT'])
    PCentry.dims[0].label='lat_bnd'
    PCentry.dims[1].label='lon_bnd'

    PCentry=f.create_dataset('CTH',data=mean['CTH'])
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

if __name__=='__main__':
    MOD03_path = '/umbc/xfs1/jianwu/users/charaj1/CMAC/MODIS-Aggregation/input-data/MYD03/'
    MOD06_path = '/umbc/xfs1/jianwu/users/charaj1/CMAC/MODIS-Aggregation/input-data/MYD06/'
    MOD06F='MYD06_L2.A2008001.0005.006.2013341193207.hdf'
    MOD03F='MYD03.A2008001.0005.006.2012066122516.hdf'
    MOD06_file=MOD06_path+MOD06F
    MOD03_file=MOD03_path+MOD03F
    lat,lon,data=read_MODIS_level2_dataV2(MOD06_file,MOD03_file)  
    CM  = data['CM']
    CTP = data['CTP']
    CTT = data['CTT']
    CTH = data['CTH']



