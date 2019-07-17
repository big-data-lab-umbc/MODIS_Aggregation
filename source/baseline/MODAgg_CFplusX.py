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

from MODAgg_IO import read_MODIS_CFplusX, save_hdfCFplusX

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
    mode = str(sys.argv[2]) 
    '''
    test: Only use 3 file couples, oneDay: 288 files
    parMonth: Parallel (a month)
    '''
    Xname = ('CTP','cloud_top_pressure_1km','hPa')
    day = int(sys.argv[3])
    import itertools
    if mode=='test':
        MOD03_path = 'input-data/MYD03/'
        MOD06_path = 'input-data/MYD06/'
    else:
        MOD03_path = '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD03/'
        MOD06_path = '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD06_L2/'

    satellite = 'Aqua'

    yr = [2008]
    mn = [1] #np.arange(1,13)  #[1]
    dy = [day] #np.arange(1,32) # [1] #np.arange(1,31)
    # latitude and longtitude boundaries of level-3 grid
    out_name=out_name+"%d%02d%02d.hdf5"%(yr[0],mn[0],day)
    lat_bnd = np.arange(-90,91,1)
    lon_bnd = np.arange(-180,180,1)
    nlat = 180
    nlon = 360

    TOT_pix      = np.zeros(nlat*nlon)
    CLD_pix      = np.zeros(nlat*nlon)
    CTP_pix      = np.zeros(nlat*nlon)

    mean = {}
    tot_F=0 #Total number of file couples read
    start=time.time()
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
#                print(MOD06_fn)
                tot_F+=1
                Lat,Lon,data = read_MODIS_CFplusX(MOD06_path+MOD06_fn,MOD03_path+MOD03_fn,Xname=Xname)
                CM  = data['CM'].ravel()
                CTP = data[Xname[0]].ravel()
                Lat=Lat.ravel()
                Lon=Lon.ravel()
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
                    CTP_pix[j] = CTP_pix[j]+np.sum(CTP[np.where(latlon_index == j)])
            granule_time += datetime.timedelta(minutes=5)

    print('derive the averaged Level-3 cloud fraction')
    print('Mode: '+mode)
    print(Xname)
    total_cloud_fraction  =  division(CLD_pix,TOT_pix).reshape([nlat,nlon])
    mean[Xname[0]] = division(CTP_pix,CLD_pix).reshape([nlat,nlon])
    save_hdfCFplusX(out_name,total_cloud_fraction,lat_bnd,lon_bnd,X=mean['CTP'],Xname=Xname)
    end=time.time()
    print('Total # of file couples read:%d'%(tot_F))
    print('Total Cloud Fraction: %0.2f'%np.nansum(total_cloud_fraction))
    print('Mean '+Xname[0]+': %0.2f '%np.nanmean(mean[Xname[0]])+Xname[2])
    print('Time elapsed:%0.2f hours'%(end/60/60-start/60/60))
   
    
#    print('plot global map')
#    plot_global_map(lat_bnd,lon_bnd,total_cloud_fraction, cmap= plt.get_cmap('jet'), \
#            vmin=0.0,vmax=1.0,title='cloud fraction', figure_name='MODIS_total_cloud_fraction_daily_mean_Python')

