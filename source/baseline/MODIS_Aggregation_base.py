#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Tue Jul 23 15:09:19 2019
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To aggregate MODIS level 2 data.
-----------------------------------------------------
value_locate(), division() from MODAgg_daily_mean.py
"""
from netCDF4 import Dataset
from jdcal import gcal2jd
import numpy as np
import time,itertools,datetime,os,sys,fnmatch

def value_locate(refx, x):
    """
    VALUE_LOCATE locates the positions of given values within a
    reference array.  The reference array need not be regularly
    spaced.  This is useful for various searching, sorting and
    interpolation algorithms.
    The reference array should be a monotonically increasing or
    decreasing list of values which partition the real numbers.  A
    reference array of NBINS numbers partitions the real number line
    into NBINS+1 regions, like so:
        REF:           X[0]         X[1]   X[2] X[3]     X[NBINS-1]
        <----------|-------------|------|---|----...---|--------------->
        INDICES:  -1           0          1    2       3        NBINS-1
        VALUE_LOCATE returns which partition each of the VALUES falls
        into, according to the figure above.  For example, a value between
        X[1] and X[2] would return a value of 1.  Values below X[0] return
        -1, and above X[NBINS-1] return NBINS-1.  Thus, besides the value
        of -1, the returned INDICES refer to the nearest reference value
        to the left of the requested value.

        Example:
            >>> refx = [2, 4, 6, 8, 10]
            >>> x = [-1, 1, 2, 3, 5, 5, 5, 8, 12, 30]
            >>> print value_locate(refx, x)
            array([-1, -1,  0,  0,  1,  1,  1,  3,  4,  4])

            This implementation is likely no the most efficient one, as there
            is
            a loop over all x, which will in practice be long. As long as x is
            shorter than 1e6 or so elements, it should still be fast (~sec).

    """

    refx = np.array(refx)
    x = np.array(x)
    loc = np.zeros(len(x), dtype='int')

    for i in np.arange(len(x)):
        ix = x[i]
        ind = ((refx - ix) <= 0).nonzero()[0]
        if len(ind) == 0:
            loc[i] = -1
        else: loc[i] = ind[-1]

    return loc

def division(n, d):

    div = np.zeros(len(d))
    for i in np.arange(len(d)):
        if d[i] >0:
          div[i]=n[i]/d[i]
        else: div[i]=None 

    return div

def readEntry(key,ncf):
    ncf.variables[key][:]
    rdval=ncf.variables[key][:]
    scale=ncf.variables[key].getncattr('scale_factor')
    offst=ncf.variables[key].getncattr('add_offset')
    return (rdval+offst)*scale

class MODIS_Level2(object):
    def __init__(self,variables):
        '''
        variables (Dictionary): {'Acronym1':('Full_name1','units1'),...}
        '''
        self.var=variables
        return
    def readHDF(self,MYD03,MYD06):
        '''
        MOD06,MOD03: Filenames including the full path.
        '''
        data={}
        myd06 = Dataset(MYD06, "r")
        CM1km = readEntry('Cloud_Mask_1km',myd06)             #Cloud mask
        data['CM'] = (np.array(CM1km[:,:,0],dtype='byte') & 0b00000110) >>1
        for key in self.var:
            data[key]=readEntry(self.var[key][0],myd06)
        
        myd03 = Dataset(MYD03, "r")
        latitude = myd03.variables["Latitude"][:,:] # Reading Specific Variable 'Latitude'.
        latitude = np.array(latitude).byteswap().newbyteorder() # Addressing Byteswap For Big Endian Error.
        longitude = myd03.variables["Longitude"][:,:] # Reading Specific Variable 'Longitude'.
        longitude = np.array(longitude).byteswap().newbyteorder() # Addressing Byteswap For Big Endian Error.
        return latitude,longitude,data

if __name__=='__main__':
    mode='test'# Only 3 files
    out_name='MODAgg_day_'+mode
    date='01/01/2008'
#        CTP = ('cloud_top_pressure_1km',myd06)     #Cloud Top Pressure (hPa)
#        CTT = ('cloud_top_temperature_1km',myd06)  #Cloud Top Temperature (K)
#        CTH = ('cloud_top_height_1km',myd06)       #Cloud Top Height (m)
    variables={'CTP':('cloud_top_pressure_1km','hPa'),'CTT':('cloud_top_temperature_1km','K')}
    MOD03_path='/home/cpnhere/cmac/input-data/MYD03/'
    MOD06_path='/home/cpnhere/cmac/input-data/MYD06/'
    
    dt=date.split('/')
    yr = [int(dt[2])]
    mn = [int(dt[0])] #np.arange(1,13)  #[1]
    dy = [int(dt[1])] #np.arange(1,32) # [1] #np.arange(1,31)
    # latitude and longtitude boundaries of level-3 grid
    out_name=out_name+"%d%02d%02d.hdf5"%(yr[0],mn[0],dy[0])
    lat_bnd = np.arange(-90,91,1)
    lon_bnd = np.arange(-180,180,1)
    nlat = 180
    nlon = 360

    TOT_pix      = np.zeros(nlat*nlon)#To compute CF 
    CLD_pix      = np.zeros(nlat*nlon)#Will be needed for all others including CF
    XXX_pix={}
    for key in variables:
        XXX_pix[key]=np.zeros(nlat*nlon)
    mean = {}
    tot_F=0 #Total number of file couples read
    start=time.time()
    MD=MODIS_Level2(variables)
    for y,m,d in  itertools.product(yr,mn, dy):
        #-------------find the MODIS prodcts--------------#
        date = datetime.datetime(y,m,d)
        JD01, JD02 = gcal2jd(y,1,1)
        JD1, JD2 = gcal2jd(y,m,d)
        JD = np.int((JD2+JD1)-(JD01+JD02) + 1)
        granule_time = datetime.datetime(y,m,d,0,0)
        while granule_time <= datetime.datetime(y,m,d,23,55):  # 23,55
#            print('granule time:',granule_time)
            MOD03_fp = 'MYD03.A{:04d}{:03d}.{:02d}{:02d}.006.?????????????.hdf'.format(y,JD,granule_time.hour,granule_time.minute)
            MOD06_fp = 'MYD06_L2.A{:04d}{:03d}.{:02d}{:02d}.006.?????????????.hdf'.format(y,JD,granule_time.hour,granule_time.minute)
            MOD03_fn, MOD06_fn =[],[]
            for MOD06_flist in  os.listdir(MOD06_path):
                if fnmatch.fnmatch(MOD06_flist, MOD06_fp):
                    MOD06_fn = MOD06_flist
            for MOD03_flist in  os.listdir(MOD03_path):
                if fnmatch.fnmatch(MOD03_flist, MOD03_fp):
                    MOD03_fn = MOD03_flist
            if MOD03_fn and MOD06_fn: # if both MOD06 and MOD03 products are in the directory
                tot_F+=1
                Lat,Lon,data = MD.readHDF(MOD03_path+MOD03_fn,MOD06_path+MOD06_fn)
                for key in data:
                    data[key]=data[key].ravel()    
                Lat=Lat.ravel()
                Lon=Lon.ravel()
                lat_index = value_locate(lat_bnd,Lat)
                lon_index = value_locate(lon_bnd,Lon)
                latlon_index = lat_index*nlon + lon_index
                latlon_index_unique = np.unique(latlon_index)
                
                for i in np.arange(latlon_index_unique.size):
                    j=latlon_index_unique[i]
                    TOT_pix[j] = TOT_pix[j]+np.sum(data['CM'][np.where(latlon_index == j)]>=0)
                    CLD_pix[j] = CLD_pix[j]+np.sum(data['CM'][np.where(latlon_index == j)]<=1) 
                    for key in data:
                        if key!='CM':
                            XXX_pix[key][j]=XXX_pix[key][j]+np.sum(data[key][np.where(latlon_index == j)])
#                    CTP_pix[j] = CTP_pix[j]+np.sum(CTP[np.where(latlon_index == j)])
#                    CTT_pix[j] = CTT_pix[j]+np.sum(CTT[np.where(latlon_index == j)])
#                    CTH_pix[j] = CTH_pix[j]+np.sum(CTH[np.where(latlon_index == j)])
            granule_time += datetime.timedelta(minutes=5)

    print('derive the averaged Level-3 cloud fraction')
    total_cloud_fraction  =  division(CLD_pix,TOT_pix).reshape([nlat,nlon])
    for key in data:
        if key!='CM':
            mean[key]=division(XXX_pix[key],CLD_pix).reshape([nlat,nlon])
