import pandas as pd
import numpy as np
import dask.array as da
import dask.dataframe as dd
import dask.delayed as delayed
import time
import math
#import graphviz
from netCDF4 import Dataset
import os,datetime,sys,fnmatch
import h5py
import dask
from dask_jobqueue import SLURMCluster
from dask.distributed import as_completed
from dask.distributed import Client
from dask.distributed import wait
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
    ncfile.close()
    return lat,lon,CM

def countzero(x, axis=1):
    #print(x)
    count0 = 0
    count1 = 0
    for i in x:
        if i <= 1:
            count0 +=1
    #print(count0/len(x))
    return (count0/len(x))

MYD06_dir= '/home/dprakas1/jianwu_common/MODIS_Aggregation/deepak-code/test/'
MYD06_prefix = 'MYD06_L2.A2008'
MYD03_dir= '/home/dprakas1/jianwu_common/MODIS_Aggregation/deepak-code/test/'
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

var_list = ['Scan Offset','Track Offset','Height Offset', 'Height', 'SensorZenith', 
            'Range', 'SolarZenith', 'SolarAzimuth', 'Land/SeaMask','WaterPresent','gflags',
            'Scan number', 'EV frames', 'Scan Type', 'EV start time', 'SD start time',
            'SV start time', 'EV center time', 'Mirror side', 'SD Sun zenith', 'SD Sun azimuth',
            'Moon Vector','orb_pos', 'orb_vel', 'T_inst2ECR', 'attitude_angles', 'sun_ref',
            'impulse_enc', 'impulse_time', 'thermal_correction', 'SensorAzimuth']



# Read Level-2 MODIS data
lat,lon,CM = read_MODIS(fname1[j],fname2[j])
if __name__ == '__main__':

    cluster = SLURMCluster(cores=32, memory='100 GB', project='pi_jianwu', queue='high_mem', walltime='02:00:00', job_extra=['--exclusive', '--qos=medium+'])

    cluster.scale(8)

    client = Client(cluster)

    dask.config.set(temporary_directory='/home/dprakas1/jianwu_common/MODIS_Aggregation/deepak-code/test-tmp')

    t0 = time.time()

    b1=[]
    cloud_pix = np.zeros((180, 360))
    for x,y in zip(fname1,fname2):
        cm = np.zeros((2030,1354), dtype=np.float32)
        lat = np.zeros((2030,1354), dtype=np.float32)
        lon = np.zeros((2030,1354), dtype=np.float32)

        #myd06 = Dataset(x, "r")
        #CM = myd06.variables["Cloud_Mask_1km"][:,:,0]# Reading Specific Variable 'Cloud_Mask_1km'.
        #CM = (np.array(CM,dtype='byte') & 0b00000110) >>1
        #CM = np.array(CM).byteswap().newbyteorder()
        d06 = xr.open_dataset(x, drop_variables="Scan Type")['Cloud_Mask_1km'][:,:,0].values
        ds06_decoded = (np.array(d06, dtype='byte') & 0b00000110) >>1
        CM = np.array(ds06_decoded).byteswap().newbyteorder()
 
        #print("CM intial shape:",CM.shape)
        cm = da.concatenate((cm,CM),axis=0)
        #print("CM shape after con:",cm.shape)
        cm=da.ravel(cm)
        #print("cm shape after ravel:",cm.shape)

        #myd03 = Dataset(y, "r")
        #latitude = myd03.variables["Latitude"][:,:]
        #longitude = myd03.variables["Longitude"][:,:]

        latitude = xr.open_dataset(y,drop_variables=var_list)['Latitude'][:,:].values
        longitude = xr.open_dataset(y,drop_variables=var_list)['Longitude'][:,:].values
        lat = da.concatenate((lat,latitude),axis=0)
        lon = da.concatenate((lon,longitude),axis=0)
        lat=da.ravel(lat)
        lon=da.ravel(lon)

        cm=cm.astype(int)
        lon=lon.astype(int)
        lat=lat.astype(int)
        lat=lat+90
        lon=lon+180
        Lat=lat.to_dask_dataframe()
        Lon=lon.to_dask_dataframe()
        CM=cm.to_dask_dataframe()
        df=dd.concat([Lat,Lon,CM],axis=1,interleave_partitions=False)


        cols = {0:'Latitude',1:'Longitude',2:'CM'}
        df = df.rename(columns=cols)

        df2=(df.groupby(['Longitude','Latitude']).CM.apply(countzero).reset_index())
        #print(type(df2))
        #df3=df2.compute()

        #df4=[df3['Longitude'].values,df3['Latitude'].values,df3['CM'].values]
        #print(df4)
        #df3=df2.compute()
        b1.append(df2)

        def results(concat_list):
            b2=dd.concat(b1)
            b_2=b2.groupby(['Latitude','Longitude']).mean().reset_index()
            combs=[]
            for x in range(0,180):
                for y in range(0,360):
                    combs.append((x, y))
            df_1=pd.DataFrame(combs)
            df_1.columns=['Latitude','Longitude']
            df_2=dd.from_pandas(df_1,npartitions=10500)
            df4=dd.merge(df_2, b_2,on=('Latitude','Longitude'), how='left')
            a = df4['CM'].to_dask_array(lengths=True)
            arr = da.from_array(a, chunks=(9257))
            final_array=arr.reshape(180,360)
            return final_array


final_arr=results(b1)

final_arr.compute()



client.close()

t1 = time.time()
total = t1-t0
print("total execution time (Seconds):" + str(total))
