import pandas as pd
import numpy as np
import dask.array as da
import dask.dataframe as dd
import dask.delayed as delayed
import time
import math
import itertools
from netCDF4 import Dataset
import os,datetime,sys,fnmatch
import matplotlib.pyplot as plt
import dask
from dask_jobqueue import SLURMCluster
from dask.distributed import as_completed
from dask.distributed import Client
from dask.distributed import wait


def read_MODIS_level2_data(MOD06_file,MOD03_file):
    #print('Reading The Cloud Mask From MOD06_L2 Product:')
    myd06 = Dataset(MOD06_file, "r")
    CM = myd06.variables["Cloud_Mask_1km"][:,:,:] # Reading Specific Variable 'Cloud_Mask_1km'.
    CM   = (np.array(CM[:,:,0],dtype='byte') & 0b00000110) >>1
    CM = np.array(CM).byteswap().newbyteorder()
   # print('The Level-2 Cloud Mask Array Shape',CM.shape)
    print(' ')

    myd03 = Dataset(MOD03_file, "r")
    #print('Reading The Latitude-Longitude From MOD03 Product:')
    latitude = myd03.variables["Latitude"][:,:] # Reading Specific Variable 'Latitude'.
    latitude = np.array(latitude).byteswap().newbyteorder() # Addressing Byteswap For Big Endian Error.
    longitude = myd03.variables["Longitude"][:,:] # Reading Specific Variable 'Longitude'.
    longitude = np.array(longitude).byteswap().newbyteorder() # Addressing Byteswap For Big Endian Error
    #print('The Level2 Latitude-Longitude Array Shape',latitude.shape)
    print(' ')

    return latitude,longitude,CM


def countzero(x, axis=1):
    #print(x)
    count0 = 0
    count1 = 0
    for i in x:
        if i <= 1:
            count0 +=1
    #print(count0/len(x))
    return count0/len(x)

MOD03_path ='/home/dprakas1/jianwu_common/MODIS_Aggregation/deepak-code/test/' #'input-data/MYD3'
MOD06_path ='/home/dprakas1/jianwu_common/MODIS_Aggregation/deepak-code/test/' #'input-data/MYD6'


MOD03_filepath = 'MYD03.A*.hdf'
MOD06_filepath = 'MYD06_L2.A*.hdf'
MOD03_filename, MOD06_filename =[],[]
MOD03_filename2, MOD06_filename2 =[],[]

for MOD06_filelist in  os.listdir(MOD06_path):
    if fnmatch.fnmatch(MOD06_filelist, MOD06_filepath):
        MOD06_filename = MOD06_filelist
        MOD06_filename2.append(MOD06_filelist)

for MOD03_filelist in  os.listdir(MOD03_path):
    if fnmatch.fnmatch(MOD03_filelist, MOD03_filepath):
        MOD03_filename = MOD03_filelist
        MOD03_filename2.append(MOD03_filelist)

#if MOD03_filename and MOD06_filename:
#    print('Reading Level 2 GeoLocation & Cloud Data')
#    Lat,Lon,CM = read_MODIS_level2_data(MOD06_path+MOD06_filename,MOD03_path+MOD03_filename)

print('The Number Of Files In The MODO3 List: ')
print(len(MOD03_filename2))
print(' ')
print('The Number Of Files In The MODO6_L2 List: ')
print(len(MOD06_filename2))
print(' ')

if __name__ == '__main__':

    cluster = SLURMCluster(cores=1, memory='100GB', project='pi_jianwu', queue='high_mem', walltime='02:00:00', job_extra=['--exclusive', '--qos=medium+'])

    cluster.scale(3)

    client = Client(cluster)

    dask.config.set(temporary_directory='/home/dprakas1/jianwu_common/MODIS_Aggregation/deepak-code/test-tmp')

    t0 = time.time()


    cm = np.zeros((2030,1354), dtype=np.float32)

    for MOD06_file in MOD06_filename2:
        MOD06_file2 = MOD06_path + MOD06_file
        myd06 = Dataset(MOD06_file2, "r")
        CM = myd06.variables["Cloud_Mask_1km"][:,:,0]# Reading Specific Variable 'Cloud_Mask_1km'.
        CM = (np.array(CM,dtype='byte') & 0b00000110) >>1
        CM = np.array(CM).byteswap().newbyteorder()
        cm = np.concatenate((cm,CM))

    print('The Cloud Mask Array Shape Is: ',cm.shape)

    lat = np.zeros((2030,1354), dtype=np.float32)
    lon = np.zeros((2030,1354), dtype=np.float32)
    for MOD03_file in MOD03_filename2:
        MOD03_file2 = MOD03_path + MOD03_file
        myd03 = Dataset(MOD03_file2, "r")
        latitude = myd03.variables["Latitude"][:,:] # Reading Specific Variable 'Latitude'.
        lat = np.concatenate((lat,latitude))
        longitude = myd03.variables["Longitude"][:,:] # Reading Specific Variable 'Longitude'.
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

    #df4=pd.merge(df_1, df3,on=('Longitude','Latitude'), how='left')

    #df5=df4['CM'].values

    #b=df5.reshape(180,360)

    client.close()

    t1 = time.time()
    total = t1-t0
    print("total execution time (Seconds):" + str(total))

    #plt.figure(figsize=(14,7))
    #plt.contourf(range(-180,180), range(-90,90), b, 100, cmap = "jet")
    #plt.xlabel("Longitude", fontsize = 14)
    #plt.ylabel("Latitude", fontsize = 14)
    #plt.title("Level 3 Cloud Fraction Aggregation for January 2008", fontsize = 16)
    #plt.colorbar()
    #plt.savefig("pandas_oneday.png")
