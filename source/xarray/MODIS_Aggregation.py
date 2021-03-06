
import dask
import numpy as np
import xarray as xr
import glob
import matplotlib.pyplot as plt
import time


class MODIS_L2_L3_Aggregation:

    def CloudFraction_Daily_Aggregation(self, satellite = "Aqua", date = "2008001"):
        var_list = ['Scan Offset','Track Offset','Height Offset', 'Height', 'SensorZenith', 'SensorAzimuth',
            'Range', 'SolarZenith', 'SolarAzimuth', 'Land/SeaMask','WaterPresent','gflags',
            'Scan number', 'EV frames', 'Scan Type', 'EV start time', 'SD start time',
            'SV start time', 'EV center time', 'Mirror side', 'SD Sun zenith', 'SD Sun azimuth',
            'Moon Vector','orb_pos', 'orb_vel', 'T_inst2ECR', 'attitude_angles', 'sun_ref',
            'impulse_enc', 'impulse_time', 'thermal_correction']


        t0 = time.time()

        if (satellite == "Aqua"):
            collection = "MYD"
        elif (satellite == "Terra"):
            collection = "MOD"
        else:
            raise ValueError("The function only accepts Aqua or Terra for satellite parameter. Your parameter value is " + satellite)

        M03_dir = "../../input-data/" + collection + "03/"
        M03_2040 = sorted(glob.glob(M03_dir + collection + "03.A" + date + "*"))

        M06_dir = "../../input-data/" + collection + "06/"
        M06_2040 = sorted(glob.glob(M06_dir + collection + "06_L2.A" + date + "*"))

        total_pix = np.zeros((180, 360))
        cloud_pix = np.zeros((180, 360))

        d06 = xr.open_mfdataset(M06_2040, concat_dim='None', parallel=True)['Cloud_Mask_1km'][:,:,:,0].values
        d06CM = d06[:,::3,::3]
        ds06_decoded = (np.array(d06CM, dtype = "byte") & 0b00000110) >> 1
        d03_lat = xr.open_mfdataset(M03_2040, concat_dim='None',drop_variables=var_list, parallel=True)['Latitude'][:,:,:].values
        d03_lon = xr.open_mfdataset(M03_2040, concat_dim='None',drop_variables=var_list, parallel=True)['Longitude'][:,:,:].values

        lat = (d03_lat[:,::3,::3].ravel() + 89.5).astype(int)
        lon = (d03_lon[:,::3,::3].ravel() + 179.5).astype(int)
        lat = np.where(lat > -1, lat, 0)
        lon = np.where(lon > -1, lon, 0)

        for i, j in zip(lat, lon):
            #print(i,j)
            total_pix[i,j] += 1

        #print(total_pix)

        index = np.nonzero(ds06_decoded.ravel() <= 0)

        cloud_lon = [lon[i] for i in index[0]]
        cloud_lat = [lat[i] for i in index[0]]

        for x, y in zip(cloud_lat, cloud_lon):
            #print(x,y)
            cloud_pix[x,y] += 1

        #print(cloud_pix)

        t1 = time.time()
        total = t1-t0
        print(total)

        total_pix[np.where(total_pix == 0)]=1.0
        cf = cloud_pix/total_pix
        print(cf)

        outputFile = collection + "08_D3.A" + date + "CloudFraction"

        plt.figure(figsize=(14,7))
        plt.contourf(range(-180,180), range(-90,90), cf, 100, cmap = "jet")
        plt.xlabel("Longitude", fontsize = 14)
        plt.ylabel("Latitude", fontsize = 14)
        plt.title("Level 3 Cloud Fraction Aggregation Result", fontsize = 16)
        plt.colorbar()
        plt.savefig("../../output-data/" + outputFile + ".png")
        cf1 = xr.DataArray(cf)
        cf1.to_netcdf("../../output-data/" + outputFile + ".nc")

        with open("../../output-data/" + outputFile + "-execution-time.txt", "a") as output:
            output.write('\n')
            output.write(str(total))

        return outputFile

aggregate = MODIS_L2_L3_Aggregation()
outputFile = aggregate.CloudFraction_Daily_Aggregation(satellite = "Aqua", date = "2008001")
print("output file is:" + outputFile)
