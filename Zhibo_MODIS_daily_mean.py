#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun May  16 17:16:55 2017

@author: zhibo
"""
from __future__ import division, print_function
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import os,datetime,sys,fnmatch
from jdcal import gcal2jd
#from plot_global_map import *
import math
import time

from cpn_MODAgg import read_MODIS_level2_dataV2, save_hdf

def read_MODIS_level2_data(MOD06_file,MOD03_file):
    print(MOD06_file)
    print(MOD03_file)
    print('reading the cloud mask from MOD06_L2 product')
    MOD06 = SD(MOD06_file, SDC.READ)
    CM1km = MOD06.select('Cloud_Mask_1km').get()
    CM   = (np.array(CM1km[:,:,0],dtype='byte') & 0b00000110) >>1
    print('level-2 cloud mask array shape',CM.shape)
    MOD06.end()

    MOD03 = SD(MOD03_file, SDC.READ)
    print('reading the lat-lon from MOD03 product')
    lat  = MOD03.select('Latitude').get()
    lon  = MOD03.select('Longitude').get()
    print('level-2 lat-lon array shape',lat.shape)
    MOD03.end()
    return lat,lon,CM

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

# beginning of the program
if __name__ == '__main__':
    out_name=sys.argv[1]
    import itertools
#    MOD03_path = '../zz_MODIS_aggregation/Shared_Sample/'#'/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD03/'
#    MOD06_path = '../zz_MODIS_aggregation/Shared_Sample/'#'/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD06_L2/'
    MOD03_path = '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD03/'
    MOD06_path = '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD06_L2/'

    satellite = 'Aqua'

    yr = [2008]
    mn = [1] #np.arange(1,13)  #[1]
    dy = [1] #np.arange(1,32) # [1] #np.arange(1,31)
    # latitude and longtitude boundaries of level-3 grid
    lat_bnd = np.arange(-90,91,1)
    lon_bnd = np.arange(-180,180,1)
    nlat = 180
    nlon = 360

    TOT_pix      = np.zeros(nlat*nlon)
    CLD_pix      = np.zeros(nlat*nlon)
    CER_pix      = np.zeros(nlat*nlon)

    tot_comT=0 #total time for computation
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
#                print('reading level 2 geolocation and cloud data')
                print(MOD06_fn)
                Lat,Lon,data = read_MODIS_level2_dataV2(MOD06_path+MOD06_fn,MOD03_path+MOD03_fn)
                CM = data['CM']
                CER=data['CER']
                start=time.time()
                Lat=Lat.ravel()
                Lon=Lon.ravel()
                CM=CM.ravel()
                CER=CER.ravel()
#                print('Total Number of pixels in this granule (cloud mask CM>=0)',np.sum(CM>=0))
#                print('Total Number of cloudy pixels (cloud mask CM<=1)',np.sum(CM<=1))
#                print('cloud fraction of this granule',np.sum(CM<=1)/np.sum(CM>=0))
#                print('projecting granule on level3 lat lon grids')
                lat_index = value_locate(lat_bnd,Lat)
                lon_index = value_locate(lon_bnd,Lon)
                latlon_index = lat_index*nlon + lon_index
#                print('computing simple level3 statistics')
                latlon_index_unique = np.unique(latlon_index)
#                print('this granule occupies',latlon_index_unique.size,'1x1 degree box')
                
                for i in np.arange(latlon_index_unique.size):
                    j=latlon_index_unique[i]
                    TOT_pix[j] = TOT_pix[j]+np.sum(CM[np.where(latlon_index == j)]>=0)
                    CLD_pix[j] = CLD_pix[j]+np.sum(CM[np.where(latlon_index == j)]<=1) 
                    CER_pix[j] = CER_pix[j]+np.sum(CER[np.where(latlon_index == j)])
                end=time.time()
                tot_comT=tot_comT+(end/60-start/60)
            granule_time += datetime.timedelta(minutes=5)

    print('derive the averaged Level-3 cloud fraction')
    total_cloud_fraction  =  division(CLD_pix,TOT_pix).reshape([nlat,nlon])
    mean_CER = division(CER_pix,CLD_pix).reshape([nlat,nlon])
    print(np.nansum(total_cloud_fraction))
    print(np.nanmean(mean_CER))
    print('Total time computational time excluding the file reading:%0.2f min'%(tot_comT))
   
    save_hdf(out_name,total_cloud_fraction,mean_CER,lat_bnd,lon_bnd)
#    print('plot global map')
#    plot_global_map(lat_bnd,lon_bnd,total_cloud_fraction, cmap= plt.get_cmap('jet'), \
#            vmin=0.0,vmax=1.0,title='cloud fraction', figure_name='MODIS_total_cloud_fraction_daily_mean_Python')

