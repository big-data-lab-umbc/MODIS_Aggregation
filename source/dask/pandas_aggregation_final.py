import pandas as pd
import numpy as np
import dask.array as da
import dask.dataframe as dd
import time
import math
from netCDF4 import Dataset
import os,datetime,sys,fnmatch
import h5py
import matplotlib.pyplot as plt
import dask
from dask_jobqueue import SLURMCluster
from dask.distributed import as_completed
from dask.distributed import Client
from dask.distributed import wait
from multiprocessing.pool import ThreadPool
import xarray as xr

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
    ncfile.close()

    ncfile=Dataset(fname2,'r')
    lat  = np.array(ncfile.variables['Latitude'])
    lon  = np.array(ncfile.variables['Longitude'])
    attr_lat = ncfile.variables['Latitude']._FillValue
    attr_lon = ncfile.variables['Longitude']._FillValue
    return lat,lon,CM

def countzero(x, axis=1):
    #print(x)
    count0 = 0
    count1 = 0
    for i in x:
        if i <= 1:
            count0 +=1
    #print(count0/len(x))
    return count0/len(x)

MYD06_dir= '/home/dprakas1/jianwu_common/MODIS_Aggregation/MODIS_one_day_data/'
MYD06_prefix = 'MYD06_L2.A2008'
MYD03_dir= '/home/dprakas1/jianwu_common/MODIS_Aggregation/MODIS_one_day_data/'

MYD03_prefix = 'MYD03.A2008'
fileformat = 'hdf'

fname1,fname2 = [],[]


days = np.arange(1,2,dtype=np.int)
for day in days:
    dc ='%03i' % day
    fname_tmp1 = read_filelist(MYD06_dir,MYD06_prefix,dc,fileformat)
    fname_tmp2 = read_filelist(MYD03_dir,MYD03_prefix,dc,fileformat)
    fname1 = np.append(fname1,fname_tmp1)
    fname2 = np.append(fname2,fname_tmp2)

# Initiate the number of day and total cloud fraction
files  = np.arange(len(fname1))



for j in range(0,1):#hdfs:
    ('steps: ',j+1,'/ ',(fname1))

    # Read Level-2 MODIS data
    #lat,lon,CM = read_MODIS(fname1[j],fname2[j])
if __name__ == '__main__':

    cluster = SLURMCluster(cores=32, memory='100GB', project='pi_jianwu', queue='high_mem', walltime='02:00:00', job_extra=['--exclusive', '--qos=medium+'])

    cluster.scale(8)

    client = Client(cluster)
    #pool= ThreadPool()
    #dask.config.set(pool=pool)
    dask.config.set(temporary_directory='/home/dprakas1/jianwu_common/MODIS_Aggregation/deepak-code/test-tmp')

    var_list = ['Scan Offset','Track Offset','Height Offset', 'Height', 'SensorZenith', 
            'Range', 'SolarZenith', 'SolarAzimuth', 'Land/SeaMask','WaterPresent','gflags',
            'Scan number', 'EV frames', 'Scan Type', 'EV start time', 'SD start time',
            'SV start time', 'EV center time', 'Mirror side', 'SD Sun zenith', 'SD Sun azimuth',
            'Moon Vector','orb_pos', 'orb_vel', 'T_inst2ECR', 'attitude_angles', 'sun_ref',
            'impulse_enc', 'impulse_time', 'thermal_correction', 'SensorAzimuth']


    t0 = time.time()

    cm = np.zeros((2030,1354), dtype=np.float32)

    for MOD06_file in fname1:
        #myd06 = Dataset(MOD06_file, "r")
        #CM = myd06.variables["Cloud_Mask_1km"][:,:,0]# Reading Specific Variable 'Cloud_Mask_1km'.
        #CM = (np.array(CM,dtype='byte') & 0b00000110) >>1
        #CM = np.array(CM).byteswap().newbyteorder()
        d06 = xr.open_dataset(MOD06_file, drop_variables="Scan Type")['Cloud_Mask_1km'][:,:,0].values
        #d06CM = d06[::3,::3]
        ds06_decoded = (np.array(d06, dtype='byte') & 0b00000110) >>1
        CM = np.array(ds06_decoded).byteswap().newbyteorder()

        cm = np.concatenate((cm,CM))

    print('The Cloud Mask Array Shape Is: ',cm.shape)

    lat = np.zeros((2030,1354), dtype=np.float32)
    lon = np.zeros((2030,1354), dtype=np.float32)
    for MOD03_file in fname2:
        #print(MOD03_file)
        #myd03 = Dataset(MOD03_file, "r")
        latitude = xr.open_dataset(MOD03_file,drop_variables=var_list)['Latitude'][:,:].values
        longitude = xr.open_dataset(MOD03_file,drop_variables=var_list)['Longitude'][:,:].values

        #latitude = myd03.variables["Latitude"][:,:] # Reading Specific Variable 'Latitude'.
        lat = np.concatenate((lat,latitude))
        #longitude = myd03.variables["Longitude"][:,:] # Reading Specific Variable 'Longitude'.
        lon = np.concatenate((lon,longitude))

    print('Longitude Shape Is: ',lon.shape)
    print('Latitude Shape Is: ',lat.shape)

    cm=np.ravel(cm)
    lat=np.ravel(lat)
    lon=np.ravel(lon)

    lon=lon.astype(int)
    lat=lat.astype(int)
    cm=cm.astype(int)
   
    d = {'Latitude' :  pd.Series(lat), 'Longitude' : pd.Series(lon),'CM':pd.Series(cm)} #Reading the values in a pandas dataframe
    df = pd.DataFrame(d,columns=['Latitude', 'Longitude','CM']) #Creating a pandas dataframe

    df2=df.groupby(['Longitude','Latitude']).CM.apply(countzero)
    df3=df2.reset_index()

    combs=[]
    for x in range(-89,91):
        for y in range(-179,181):
            combs.append((x, y))

    df_1=pd.DataFrame(combs)
    df_1.columns=['Latitude','Longitude']

    df4=pd.merge(df_1, df3,on=('Longitude','Latitude'), how='left')

    df5=df4['CM'].values

    b=df5.reshape(180,360)

    client.close()

    t1 = time.time()
    total = t1-t0
    print("total execution time (Seconds):" + str(total))

    plt.figure(figsize=(14,7))
    plt.contourf(range(-180,180), range(-90,90),b, 100, cmap = "jet")
    plt.xlabel("Longitude", fontsize = 14)
    plt.ylabel("Latitude", fontsize = 14)
    plt.title("Level 3 Cloud Fraction Day level Aggregation for Dask", fontsize = 16)
    plt.colorbar()
    plt.savefig("dask_aggregation_day_level.png")
