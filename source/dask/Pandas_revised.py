#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import dask.array as da
import dask.dataframe as dd
import time
import math
#import graphviz
from netCDF4 import Dataset
import os,datetime,sys,fnmatch
import h5py


def read_filelist(loc_dir,prefix,unie,fileformat):
    # Read the filelist in the specific directory
    str = os.popen("ls "+ loc_dir + prefix + unie + "*."+fileformat).read()
    fname = np.array(str.split("\n"))
    fname = np.delete(fname,len(fname)-1)
    
    return fname

def read_MODIS(fname1,fname2,verbose=False): # READ THE HDF FILE
    # Read the cloud mask from MYD06_L2 product')
    ncfile=Dataset(fname1,'r')
    CM1km = np.array(ncfile.variables['Cloud_Mask_1km'])
    CM   = (np.array(CM1km[:,:,0],dtype='byte') & 0b00000110) >>1
    #ncfile = Dataset(fname1, "r")
    #CM = myd06.variables["Cloud_Mask_1km"][:,:,:] # Reading Specific Variable 'Cloud_Mask_1km'.
    #CM   = (np.array(CM[:,:,0],dtype='byte') & 0b00000110) >>1
    ncfile.close()
    
    ncfile=Dataset(fname2,'r')
    lat  = np.array(ncfile.variables['Latitude'])
    lon  = np.array(ncfile.variables['Longitude'])
    #ncfile = Dataset(MOD03_files, "r")
    #latitude = myd03.variables["Latitude"][:,:] # Reading Specific Variable 'Latitude'.
    #latitude = np.array(latitude).byteswap().newbyteorder() # Addressing Byteswap For Big Endian Error.
    #longitude = myd03.variables["Longitude"][:,:] # Reading Specific Variable 'Longitude'.
    attr_lat = ncfile.variables['Latitude']._FillValue
    attr_lon = ncfile.variables['Longitude']._FillValue
    """#Use _FillValue to remove fill data in lat & lon
    lat[np.where(lat == attr_lat)] = 0.0
    lon[np.where(lat == attr_lat)] = 0.0
    CM [np.where(lat == attr_lat)] = 0.5 #which will not be identified by lines 80-83 
    lat[np.where(lon == attr_lon)] = 0.0
    lon[np.where(lon == attr_lon)] = 0.0
    CM [np.where(lon == attr_lon)] = 0.5 #which will not be identified by lines 80-83
    ncfile.close()"""
    return lat,lon,CM

#Function for processing the cloud Fraction
def value_locate(refx, x):
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

satellite = 'Aqua'

yr = [2008]
mn = [1] #np.arange(1,13)  #[1]
dy = [1] #np.arange(1,32) # [1] #np.arange(1,31)
lat_bnd = np.arange(-90,91,1)# latitude and longtitude boundaries of level-3 grid
lon_bnd = np.arange(-180,180,1)
nlat = 180
nlon = 360

TOT_pix      = np.zeros(nlat*nlon)
CLD_pix      = np.zeros(nlat*nlon)

MYD06_dir= '/home/dprakas1/jianwu_common/MODIS_Aggregation/MODIS_one_day_data/'
MYD06_prefix = 'MYD06_L2.A2008'
MYD03_dir= '/home/dprakas1/jianwu_common/MODIS_Aggregation/MODIS_one_day_data/'
MYD03_prefix = 'MYD03.A2008'
fileformat = 'hdf'

#MYD06_dir= '/home/dprakas1/jianwu_common/MODIS_Aggregation/MODIS_one_day_data-2/'
#MYD06_prefix = 'MYD06_L2.A2016'
#MYD03_dir= '/home/dprakas1/jianwu_common/MODIS_Aggregation/MODIS_one_day_data-2/'
#MYD03_prefix = 'MYD03.A2016'
#fileformat = 'hdf'

fname1,fname2 = [],[]


days = np.arange(1,31,dtype=np.int)
for day in days:
    dc ='%03i' % day
    fname_tmp1 = read_filelist(MYD06_dir,MYD06_prefix,dc,fileformat)
    fname_tmp2 = read_filelist(MYD03_dir,MYD03_prefix,dc,fileformat)
    fname1 = np.append(fname1,fname_tmp1)
    fname2 = np.append(fname2,fname_tmp2)

# Initiate the number of day and total cloud fraction
files  = np.arange(len(fname1))



for j in range(0,1):#hdfs:
    #print('steps: ',j+1,'/ ',(fname1)) 

    # Read Level-2 MODIS data
    lat,lon,CM = read_MODIS(fname1[j],fname2[j])
#print((fname1))
#print((fname2))
#rint(CM)
#lat = lat.ravel()
#lon = lon.ravel()
#CM  = CM.ravel()
CM.shape    

cm = np.zeros((2030,1354), dtype=np.float32)
bit0 = np.zeros((2030,1354), dtype=np.float32)
bit12 = np.zeros((2030,1354), dtype=np.float32)    

for MOD06_file in fname1:
    #print(MOD06_file)
    myd06 = Dataset(MOD06_file, "r")
    CM = myd06.variables["Cloud_Mask_1km"][:,:,0]# Reading Specific Variable 'Cloud_Mask_1km'.
    CM = (np.array(CM,dtype='byte') & 0b00000110) >>1
    CM = np.array(CM).byteswap().newbyteorder()
    #print(CM.shape)
    cm = np.concatenate((cm,CM))
    #bit0 = np.dstack((bit0,bit0r))
    #bit12 = np.dstack((bit12,bit12r))
    
print('The Cloud Mask Array Shape Is: ',cm.shape)

#myd03_name = '/home/dprakas1/jianwu_common/MODIS_Aggregation/MODIS_one_day_data/'
lat = np.zeros((2030,1354), dtype=np.float32)
lon = np.zeros((2030,1354), dtype=np.float32)
for MOD03_file in fname2:
    #print(MOD03_file)
    myd03 = Dataset(MOD03_file, "r")
    latitude = myd03.variables["Latitude"][:,:] # Reading Specific Variable 'Latitude'.
    lat = np.concatenate((lat,latitude))


    longitude = myd03.variables["Longitude"][:,:] # Reading Specific Variable 'Longitude'.
    lon = np.concatenate((lon,longitude))
    
print('Longitude Shape Is: ',lon.shape)
print('Latitude Shape Is: ',lat.shape)

cm=np.ravel(cm)
lat=np.ravel(lat)
lon=np.ravel(lon)

lat = lat.astype(int)
lon = lon.astype(int)
cm=cm.astype(int)

d = {'Latitude' :  pd.Series(lat), 'Longitude' : pd.Series(lon),'CM':pd.Series(cm)} #Reading the values in a pandas dataframe
df = pd.DataFrame(d,columns=['Latitude', 'Longitude','CM']) #Creating a pandas dataframe
#df=df.astype(np.int8)
df2=df.groupby(['Longitude','Latitude']).CM.apply(countzero)

df3=df2.reset_index()
df3

import pandas as pd
combs=[]
for x in range(-89,91):
    for y in range(-179,181):
        combs.append((x, y))


df_1=pd.DataFrame(combs)
df_1.columns=['Latitude','Longitude']


df4=pd.merge(df_1, df3,on=('Longitude','Latitude'), how='left')

df5=df4['CM'].values

b=df5.reshape(180,360)


total_cloud_fraction = b
out_name = '/home/dprakas1/jianwu_common/MODIS_Aggregation/deepak-code/reshape_full_2008.hdf5'
def save_hdf(out_name,total_cloud_fraction,lat_bnd,lon_bnd):
    f=h5py.File(out_name,'w')
    PCentry=f.create_dataset('CF',data=total_cloud_fraction)
    PCentry.dims[0].label='lat_bnd'
    #PCentry.dims[1].label='lon_bnd'

    PC=f.create_dataset('lat_bnd',data=lat_bnd)
    PC.attrs['units']='degrees'
    PC.attrs['long_name']='Latitude_boundaries'

    PC=f.create_dataset('lon_bnd',data=lon_bnd)
    PC.attrs['units']='degrees'
    PC.attrs['long_name']='Longitude_boundaries'
    f.close()
    print(out_name+' Saved!!')

save_hdf(out_name,total_cloud_fraction,lat_bnd,lon_bnd)



#for one day's data
import h5py
import numpy as np
from comparisons import readData, doPlot
benchmark_p="/umbc/xfs1/jianwu/users/charaj1/CMAC/MODIS-Aggregation/output-data/benchmark/MODAgg_3var_parMonth/"
CF_BMK,_,_=readData(benchmark_p+"MODAgg_3var_parMonth_20080101.hdf5")#Benchmark

f=h5py.File('/home/dprakas1/jianwu_common/MODIS_Aggregation/deepak-code/reshape_full_2008.hdf5','r')
CF_res=f['CF'][:]
fig1,fig1_ttl=doPlot(CF_res,CF_BMK,'Comparison')

fig1.savefig('/home/dprakas1/jianwu_common/MODIS_Aggregation/deepak-code/Reshape_full_2008')
