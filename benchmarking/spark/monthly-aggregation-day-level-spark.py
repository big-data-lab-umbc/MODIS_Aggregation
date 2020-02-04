from netCDF4 import Dataset
import numpy as np
import glob
import matplotlib.pyplot as plt
import time
import h5py
from pyspark.sql import SparkSession

def aggregateOneDayData(z):

    var_list = ['Scan Offset','Track Offset','Height Offset', 'Height', 'SensorZenith', 
            'Range', 'SolarZenith', 'SolarAzimuth', 'Land/SeaMask','WaterPresent','gflags',
            'Scan number', 'EV frames', 'Scan Type', 'EV start time', 'SD start time',
            'SV start time', 'EV center time', 'Mirror side', 'SD Sun zenith', 'SD Sun azimuth',
            'Moon Vector','orb_pos', 'orb_vel', 'T_inst2ECR', 'attitude_angles', 'sun_ref',
            'impulse_enc', 'impulse_time', 'thermal_correction', 'SensorAzimuth']
            
            
    total_pix = np.zeros((180, 360))
    cloud_pix = np.zeros((180, 360))

    M03_files = sorted(glob.glob(M03_dir + "MYD03.A2008" + "z" + "*"))
    M06_files = sorted(glob.glob(M06_dir + "MYD06_L2.A2008" + "z" + "*"))

    for x,y in zip(M06_files,M03_files):

        d06 = xr.open_dataset(x, drop_variables="Scan Type")['Cloud_Mask_1km'][:, :, 0].values
        d06CM = d06[::3, ::3]
        ds06_decoded = (np.array(d06CM, dtype="byte") & 0b00000110) >> 1
        d03_lat = xr.open_dataset(y, drop_variables=var_list)['Latitude'][:, :].values
        d03_lon = xr.open_dataset(y, drop_variables=var_list)['Longitude'][:, :].values

        lat = (d03_lat[::3, ::3].ravel() + 89.5).astype(int)
        lon = (d03_lon[::3, ::3].ravel() + 179.5).astype(int)
        lat = np.where(lat > -1, lat, 0)
        lon = np.where(lon > -1, lon, 0)

        for i, j in zip(lat, lon):
            total_pix[i, j] += 1

        index = np.nonzero(ds06_decoded.ravel() == 0)

        cloud_lon = [lon[i] for i in index[0]]
        cloud_lat = [lat[i] for i in index[0]]

        for x, y in zip(cloud_lat, cloud_lon):
            cloud_pix[x, y] += 1

    return cloud_pix, total_pix

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
    
def save_output():
    cf1 = xr.DataArray(cf)
    cf1.to_netcdf("monthlyCloudFraction-day-level-parallelization.nc")
    plt.figure(figsize=(14, 7))
    plt.contourf(range(-180, 180), range(-90, 90), cf, 100, cmap="jet")
    plt.xlabel("Longitude", fontsize=14)
    plt.ylabel("Latitude", fontsize=14)
    plt.title("Level 3 Cloud Fraction Aggregation for January 2008", fontsize=16)
    plt.colorbar()
    plt.savefig("monthlyCloudFraction-day-level-parallelization.png")


if __name__ == '__main__':

    M06_dir = "/Users/jianwu/Documents/github/MODIS-Aggregation/input-data/MYD06/"
    M03_dir = "/Users/jianwu/Documents/github/MODIS-Aggregation/input-data/MYD03/"
    M06_files = sorted(glob.glob(M06_dir + "MYD06_L2.A2008*"))
    M03_files = sorted(glob.glob(M03_dir + "MYD03.A2008*"))
    file_pairs = zip(M06_files, M03_files)
    print(file_pairs)

    index = 31
    y = [str(x).zfill(3) for x in range(index + 1)]
    z = y[1:]
    print(z)

    t0 = time.time()
    
    # Initiate and process the parallel by Spark
    spark = SparkSession\
            .builder\
            .appName("MODIS_agg")\
            .getOrCreate()
    sc = spark.sparkContext
    global_cloud_pix, global_total_pix = sc.parallelize(z, 31).map(lambda x: aggregateOneDayData(x)).reduce(lambda x, y: (x[0] + y[0], x[1] + y[1]))
    spark.stop() # Stop Spark
    lat_bnd = np.arange(-90,90,1)
    lon_bnd = np.arange(-180,180,1)
    global_total_pix[np.where(global_total_pix == 0)]=1.0
    total_cloud_fraction = (global_cloud_pix/global_total_pix)
    print("total_cloud_fraction:" + str(total_cloud_fraction))
    print("total_cloud_fraction.shape:" + str(total_cloud_fraction.shape))
    #total_cloud_fraction = (global_cloud_pix/global_total_pix).reshape([lat_bnd,lon_bnd])
    save_output()

    #calculate execution time
    t1 = time.time()
    total = t1-t0
    print("total execution time (Seconds):" + str(total))
