from netCDF4 import Dataset
import numpy as np
import glob
import matplotlib.pyplot as plt
import time
import h5py
import xarray as xr
from pyspark.sql import SparkSession

def aggregateOneFileData(M06_file, M03_file):
    """Aggregate one file from MYD06_L2 and its corresponding file from MYD03. Read 'Cloud_Mask_1km' variable from the MYD06_L2 file, read 'Latitude' and 'Longitude' variables from the MYD03 file. Group Cloud_Mask_1km values based on their (lat, lon) grid.
    Args:
        M06_file (string): File path for M06_file.
        M03_file (string): File path for corresponding M03_file.
        
    Returns:
        (cloud_pix, total_pix) (tuple): cloud_pix is an 2D(180*360) numpy array for cloud pixel count of each grid, total_pix is an 2D(180*360) numpy array for total pixel count of each grid.
    """
    #print(M06_file)
    #print(M03_file)
    
    #output as a dictionary
    output_dictionaray = {} 
    output_list = []
    total_pix = np.zeros((180, 360))
    cloud_pix = np.zeros((180, 360))
    #read 'Cloud_Mask_1km' variable from the MYD06_L2 file, whose shape is (2030, 1354)
    ncfile = Dataset(M06_file,'r')
    d06 = ncfile.variables['Cloud_Mask_1km'][:,:,0]
    ncfile.close()
    # sampling data with 1/3 ratio (pick 1st, 4th, 7th, ...) in both latitude and longitude direction. d06CM's shape is (677, 452)
    d06CM = d06[::3,::3]
    ds06_decoded = (np.array(d06CM, dtype = "byte") & 0b00000110) >> 1
    ds06_1d = ds06_decoded.ravel()
    
    
    #[print("one element:" + str(a)) for a in ds06_1d]
    # shape of d03_lat and d03_lon: (2030, 1354)
    ncfile = Dataset(M03_file,'r')
    d03_lat = ncfile.variables['Latitude'][:,:]
    d03_lon = ncfile.variables['Longitude'][:,:]
    ncfile.close()
    
    # sampling data with 1/3 ratio, shape of lat and lon: (677, 452), then convert data from 2D to 1D, then add offset to change value range from (-90, 90) to (0, 180) for lat.
    lat = (d03_lat[::3,::3].ravel() + 89.5).astype(int)
    lon = (d03_lon[::3,::3].ravel() + 179.5).astype(int)
    lat = np.where(lat > -1, lat, 0)
    lon = np.where(lon > -1, lon, 0)
    
    #create one element <(lat, lon), (cloud_pixel_number, total_pixel_number)> for each pixel and add it to output list
    for i in range(0, ds06_1d.size):
        #print("one element:" + str(lat[i]) + "," + str(lon[i]) + ", " + str(ds06_1d[i]))
        output_list.append(((lat[i], lon[i]), (1 if ds06_1d[i] == 0 else 0, 1)))
    
    #print("output for " + str(M06_file) + ":" + str(output_list))    
    return output_list


def save_output(cf):
    cf1 = xr.DataArray(cf)
    cf1.to_netcdf("daily-CloudFraction-day-level-parallelization.nc")
    plt.figure(figsize=(14, 7))
    plt.contourf(range(-180, 180), range(-90, 90), cf, 100, cmap="jet")
    plt.xlabel("Longitude", fontsize=14)
    plt.ylabel("Latitude", fontsize=14)
    plt.title("Level 3 Cloud Fraction Aggregation for January 1st 2008", fontsize=16)
    plt.colorbar()
    plt.savefig("daily-CloudFraction-grid-level-parallelization.png")

if __name__ =='__main__':

    M06_dir = "/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD06_L2/"
    M03_dir = "/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD03/"
    M06_files = sorted(glob.glob(M06_dir + "MYD06_L2.A2008*"))
    file_num = len(M06_files)
    M03_files = sorted(glob.glob(M03_dir + "MYD03.A2008*"))
    file_pairs = zip(M06_files, M03_files)
    #print(list(file_pairs))
    #print(len(list(file_pairs)))    

    t0 = time.time()
    
    # Initiate and process the parallel by Spark
    spark = SparkSession\
            .builder\
            .appName("MODIS_agg")\
            .getOrCreate()
    sc = spark.sparkContext
    result = sc.parallelize(list(file_pairs), file_num).flatMap(lambda x: aggregateOneFileData(x[0],x[1])).reduceByKey(lambda x, y: (x[0] + y[0], x[1] + y[1])).collect()
    spark.stop() # Stop Spark
    
    #print(result)
    
    global_cloud_pix = np.zeros((180, 360))
    global_total_pix = np.zeros((180, 360))
    for element in result:
        #print("element:" + str(element))
        global_cloud_pix[element[0][0], element[0][1]] = element[1][0]
        global_total_pix[element[0][0], element[0][1]] = element[1][1]
    
    total_cloud_fraction = (global_cloud_pix/global_total_pix)
    print("total_cloud_fraction:" + str(total_cloud_fraction))
    print("total_cloud_fraction.shape:" + str(total_cloud_fraction.shape))
    
    #total_cloud_fraction = (global_cloud_pix/global_total_pix).reshape([lat_bnd,lon_bnd])
    save_output(total_cloud_fraction)

    #calculate execution time
    t1 = time.time()
    total = t1-t0
    print("total execution time (Seconds):" + str(total))
