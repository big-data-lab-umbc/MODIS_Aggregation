#!/usr/bin/env python
# coding:utf8
# -*- coding: utf-8 -*-
"""
Main Program: PLOT MODIS AGGREGATION RESULT 

Created on 2019

@author: Jianyu Zheng
"""
import os 
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# files store output
#cld_fraction = np.loadtxt('cld_fraction_mpi_files_test.dat')

filename1 = "MYD08_M3A200801_baseline.h5"            #"MOD08_M3A200801_test04.h5"
filename2 = "MYD08_M3.A2008001.006.2014263192806.hdf"#"MYD08_D3.A2008001.006.2015088011904.hdf"

hdf = Dataset(filename1,'r')
cld_fraction1 = np.array(hdf.variables["Cloud_Fraction_Mean"][:])
#fillvalue     = hdf.variables["Cloud_Fraction_Mean"]._FillValue
#cld_fraction1[np.where(cld_fraction1 == fillvalue)] = 0.0

hdf = Dataset(filename2,'r')
cld_fraction2 = np.flip(hdf.variables["Cloud_Fraction_Mean_Mean"][:],0)
scale_factor  = hdf.variables["Cloud_Fraction_Mean_Mean"].scale_factor
fillvalue     = hdf.variables["Cloud_Fraction_Mean_Mean"]._FillValue
cld_fraction2[np.where(cld_fraction2 == fillvalue)] = 0.0
#cld_fraction2 = cld_fraction2*scale_factor

diff = cld_fraction1.round(4) - cld_fraction2.round(4)
#diff[np.where(((diff/cld_fraction2)*100) > 100)] = np.nan
diff[np.where(((diff/cld_fraction2)*100) < 165)] = np.nan
print(cld_fraction1,cld_fraction2,scale_factor)

lon = np.arange(-180,180,1)
lat = np.arange(-90,90,1)
Lon,Lat = np.meshgrid(lon,lat)
levels = 100

plt.figure(figsize=(10,5)) 
m = Basemap(projection='cyl',fix_aspect=False,lon_0=0,lat_0=0)
m.drawcountries()
m.drawcoastlines(linewidth=0.4)
#((diff/cld_fraction2)*100).astype(int),
cset = plt.contourf(Lon, Lat, cld_fraction2, levels, cmap='jet')#, vmin=165, vmax=100)
ax = plt.gca()
ax.set_xticks([-150,-75,0,75,150])
ax.set_yticks([-60,-30,0,30,60])
plt.xlabel("Longitude (degree)",fontsize=12)
plt.ylabel("Latitude (degree)",fontsize=12)
#plt.title('MODIS Monthly Mean Cloud Fraction Difference \n (Product - Origin)/Origin (01/2008)',fontsize=15)
plt.title('MODIS Monthly Mean Cloud Fraction from Original Product (01/2008)',fontsize=15)

cg = plt.colorbar(cset)
#cg.set_ticks([0.0,0.2,0.4,0.6,0.8,1.0]) #np.arange(165,250,6)
#cg.set_label('Cloud Total Fraction Difference (%)',fontsize=12)
cg.set_label('Cloud Total Fraction',fontsize=12)

plt.savefig('MODIS_Mean_monthly_origin.png',dpi=600)
plt.show()
plt.close()