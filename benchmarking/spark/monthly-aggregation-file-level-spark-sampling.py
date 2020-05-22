from netCDF4 import Dataset
import numpy as np
import glob
import matplotlib.pyplot as plt
import time
import h5py
import xarray as xr
import sys
from pyspark.sql import SparkSession

def aggregateOneFileData(M06_file, M03_file, sampling_rate):
    """Aggregate one file from MYD06_L2 and its corresponding file from MYD03. Read 'Cloud_Mask_1km' variable from the MYD06_L2 file, read 'Latitude' and 'Longitude' variables from the MYD03 file. Group Cloud_Mask_1km values based on their (lat, lon) grid.
    Args:
        M06_file (string): File path for M06_file.
        M03_file (string): File path for corresponding M03_file.
        
    Returns:
        (cloud_pix, total_pix) (tuple): cloud_pix is an 2D(180*360) numpy array for cloud pixel count of each grid, total_pix is an 2D(180*360) numpy array for total pixel count of each grid.
    """
    #sampling_rate = 1 # 1 meaning no sampling,
    print(M06_file)
    print(M03_file) 
    total_pix = np.zeros((180, 360))
    cloud_pix = np.zeros((180, 360))
    #read 'Cloud_Mask_1km' variable from the MYD06_L2 file, whose shape is (2030, 1354)
    d06 = Dataset(M06_file, "r").variables['Cloud_Mask_1km'][:,:,0]
    # sampling data with 1/3 ratio (pick 1st, 4th, 7th, ...) in both latitude and longitude direction. d06CM's shape is (677, 452)
    d06CM = d06[::sampling_rate,::sampling_rate]
    ds06_decoded = (np.array(d06CM, dtype = "byte") & 0b00000110) >> 1
    # shape of d03_lat and d03_lon: (2030, 1354)
    d03_lat = Dataset(M03_file, "r").variables['Latitude'][:,:]
    d03_lon = Dataset(M03_file, "r").variables['Longitude'][:,:]
    
    # sampling data with 1/3 ratio, shape of lat and lon: (677, 452), then convert data from 2D to 1D, then add offset to change value range from (-90, 90) to (0, 180) for lat.
    lat = (d03_lat[::sampling_rate,::sampling_rate].ravel() + 89.5).astype(int)
    lon = (d03_lon[::sampling_rate,::sampling_rate].ravel() + 179.5).astype(int)
    lat = np.where(lat > -1, lat, 0)
    lon = np.where(lon > -1, lon, 0)
    # increment total_pix by 1 for the grid for each value in (lat, lon).
    for i, j in zip(lat, lon):
            total_pix[i,j] += 1
    
    # covert ds06_decoded from 2D to 1D, check whether each element is less than or equal to 0, return a tuple whose first element is an 1D arrays of indices of ds06_decoded's elements whose value is less than or equal to 0.  
    index = np.nonzero(ds06_decoded.ravel() == 0)
    # get its lat and lon for each cloud pixel.
    # we can use this approach because the internal structure (677, 452) is the same for both MYD03 and MYD06.
    cloud_lon = [lon[i] for i in index[0]]
    cloud_lat = [lat[i] for i in index[0]]
     # increment cloud_pix by 1 for the grid for each value in (cloud_lat, cloud_lon).
    for x, y in zip(cloud_lat, cloud_lon):
        cloud_pix[x,y] += 1  
        
    return cloud_pix, total_pix


def save_output(cf, node_num, sampling_rate):
    cf1 = xr.DataArray(cf)
    output_file_name = "monthlyCloudFraction-file-level-" + node_num + "-nodes-" + sampling_rate + "-sampling"
    cf1.to_netcdf(output_file_name + ".nc")
    plt.figure(figsize=(14, 7))
    plt.contourf(range(-180, 180), range(-90, 90), cf, 100, cmap="jet")
    plt.xlabel("Longitude", fontsize=14)
    plt.ylabel("Latitude", fontsize=14)
    plt.title("Level 3 Cloud Fraction Aggregation for January 2008", fontsize=16)
    plt.colorbar()
    plt.savefig(output_file_name + ".png")

if __name__ =='__main__':

    #node_num = int(sys.argv[1])
    if len(sys.argv) == 2:
        sampling_rate = int(sys.argv[2])
    else:
        sampling_rate = 3 # Default sampling rate is 3
    print ("running on " + sys.argv[1] + "nodes with " + sys.argv[2] + " sampling.")
    t0 = time.time()

    #M06_dir = "/umbc/xfs1/jianwu/common/MODIS_Aggregation/MODIS_one_day_data/"
    #M03_dir = "/umbc/xfs1/jianwu/common/MODIS_Aggregation/MODIS_one_day_data/"
    M06_dir = "/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD06_L2/"
    M03_dir = "/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD03/"
    M06_files = sorted(glob.glob(M06_dir + "MYD06_L2.A2008*"))
    file_num = len(M06_files)
    M03_files = sorted(glob.glob(M03_dir + "MYD03.A2008*"))
    file_pairs = zip(M06_files, M03_files)
    #print(list(file_pairs))
    #print(len(list(file_pairs)))    

    # Initiate and process the parallel by Spark
    spark = SparkSession\
            .builder\
            .appName("MODIS_agg")\
            .getOrCreate()
    sc = spark.sparkContext
    global_cloud_pix, global_total_pix = sc.parallelize(list(file_pairs), file_num)\
                                           .map(lambda x: aggregateOneFileData(x[0], x[1], sampling_rate))\
                                           .reduce(lambda x, y: (x[0] + y[0], x[1] + y[1]))
    spark.stop() # Stop Spark
    lat_bnd = np.arange(-90,90,1)
    lon_bnd = np.arange(-180,180,1)
    global_total_pix[np.where(global_total_pix == 0)]=1.0
    total_cloud_fraction = (global_cloud_pix/global_total_pix)
    print("total_cloud_fraction:" + str(total_cloud_fraction))
    print("total_cloud_fraction.shape:" + str(total_cloud_fraction.shape))

    #calculate execution time
    t1 = time.time()
    total = t1-t0
    print("total execution time (Seconds):" + str(total))    
    
    #total_cloud_fraction = (global_cloud_pix/global_total_pix).reshape([lat_bnd,lon_bnd])
    save_output(total_cloud_fraction, sys.argv[1], sys.argv[2])
