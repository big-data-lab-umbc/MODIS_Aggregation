import pandas as pd
import numpy as np
import dask.array as da
import dask.dataframe as dd
import time
import math
from netCDF4 import Dataset
import os,datetime,sys,fnmatch
import h5py
import dask.delayed as delayed
import xarray as xr
import matplotlib.pyplot as plt
import dask
from dask_jobqueue import SLURMCluster
from dask.distributed import as_completed
from dask.distributed import Client
import xarray as xr
from dask.distributed import wait
import gc

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
    lat  = (np.array(ncfile.variables['Latitude']))
    lon  = (np.array(ncfile.variables['Longitude']))
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

if __name__ == '__main__':
    MYD06_dir= '/home/dprakas1/jianwu_common/MODIS_Aggregation/deepak-code/test/'
    MYD06_prefix = 'MYD06_L2.A2008'
    MYD03_dir= '/home/dprakas1/jianwu_common/MODIS_Aggregation/deepak-code/test/'
    MYD03_prefix = 'MYD03.A2008'

    fileformat = 'hdf'#

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

    var_list = ['Scan Offset','Track Offset','Height Offset', 'Height', 'SensorZenith',
            'Range', 'SolarZenith', 'SolarAzimuth', 'Land/SeaMask','WaterPresent','gflags',
            'Scan number', 'EV frames', 'Scan Type', 'EV start time', 'SD start time',
            'SV start time', 'EV center time', 'Mirror side', 'SD Sun zenith', 'SD Sun azimuth',
            'Moon Vector','orb_pos', 'orb_vel', 'T_inst2ECR', 'attitude_angles', 'sun_ref',
            'impulse_enc', 'impulse_time', 'thermal_correction', 'SensorAzimuth']



    cluster = SLURMCluster(cores=32, memory='100 GB',processes=True,threads_per_worker=1, project='pi_jianwu', queue='high_mem',    walltime='02:00:00', job_extra=['--exclusive', '--qos=medium+'])

    cluster.scale(8)

    client = Client(cluster)

    print(cluster.job_script())

    dask.config.set(temporary_directory='/umbc/xfs1/jianwu/common/MODIS_Aggregation/deepak-code/test-tmp-jianwu')
    pd.options.mode.chained_assignment = None

    t0 = time.time()


    cm = np.zeros((2030,1354), dtype=np.float32)

    for MOD06_file in fname1:
        #print(MOD06_file)
        #myd06 = Dataset(MOD06_file, "r")
        #CM = myd06.variables["Cloud_Mask_1km"][:,:,0]# Reading Specific Variable 'Cloud_Mask_1km'.
        #CM = (np.array(CM,dtype='byte') & 0b00000110) >>1
        #CM = (np.array(CM).byteswap().newbyteorder())
        d06 = xr.open_dataset(MOD06_file, drop_variables="Scan Type")['Cloud_Mask_1km'][:,:,0].values
        ds06_decoded = (np.array(d06, dtype='byte') & 0b00000110) >>1
        CM = np.array(ds06_decoded).byteswap().newbyteorder()


        cm = da.concatenate((cm,CM),axis=0)



        #cm = da.from_array(CM, chunks =(2030,1354))
        #cm = delayed(np.concatenate((cm,CM)))
        cm = da.concatenate((cm,CM),axis=0)

    print('The Cloud Mask Array Shape Is: ',cm.shape)


    lat = np.zeros((2030,1354), dtype=np.float32)
    lon = np.zeros((2030,1354), dtype=np.float32)

    for MOD03_file in fname2:
        #print(MOD03_file)
        #myd03 = Dataset(MOD03_file, "r")
        #latitude = myd03.variables["Latitude"][:,:]# Reading Specific Variable 'Latitude'.
        #lat = da.from_array(latitude, chunks =(2030,1354))
        latitude = xr.open_dataset(MOD03_file,drop_variables=var_list)['Latitude'][:,:].values


        lat = da.concatenate((lat,latitude),axis=0)


        #longitude = myd03.variables["Longitude"][:,:] # Reading Specific Variable 'Longitude'.
        longitude = xr.open_dataset(MOD03_file,drop_variables=var_list)['Longitude'][:,:].values
        #lon = da.from_array(longitude, chunks =(2030,1354))
        lon = da.concatenate((lon,longitude),axis=0)

    print('Longitude Shape Is: ',lon.shape)
    print('Latitude Shape Is: ',lat.shape)

    cm=da.ravel(cm)
    lat=da.ravel(lat)
    lon=da.ravel(lon)

    lon=lon.astype(int)
    lat=lat.astype(int)
    cm=cm.astype(int)

    Lat=lat.to_dask_dataframe()
    Lon=lon.to_dask_dataframe()
    CM=cm.to_dask_dataframe()

    df=dd.concat([Lat,Lon,CM],axis=1,interleave_partitions=False)

    cols = {0:'Latitude',1:'Longitude',2:'CM'}
    df = df.rename(columns=cols)
    
    df=client.persist(df)
    df2=delayed(df.groupby(['Longitude','Latitude']).CM.apply(countzero).reset_index())
    df3=df2.compute()
    #print(gc.collect())
    #print(df2)
    #df3=client.compute(df2)
    #print(df3)
    #tt=client.gather(df3)
    #print(tt)
    client.close()

    combs=[]
    for x in range(-89,91):
        for y in range(-179,181):
            combs.append((x, y))
        
    df_1=pd.DataFrame(combs)
    df_1.columns=['Latitude','Longitude']
    #df_2=dd.from_pandas(df_1,npartitions=1)

    df4=pd.merge(df_1, df3,on=('Longitude','Latitude'), how='left')
    
    
    df5=df4['CM'].values
    b=df5.reshape(180,360)
    client.close()


    t1 = time.time()
    total = t1-t0
    print("total execution time (Seconds):" + str(total))
