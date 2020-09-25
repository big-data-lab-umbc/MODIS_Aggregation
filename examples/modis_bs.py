import sys
import numpy as np
import pandas as pd
import timeit
from datetime import date, datetime
from dateutil.rrule import rrule, DAILY, MONTHLY
import h5py
from MODIS import *

if __name__ == '__main__':
    # This is the main program for using concurrent to speed up the whole process

    # -------------STEP 0: Read the input from User --------
    # checking user input
    if (len(sys.argv) != 16) & (len(sys.argv) != 17):
        print("Wrong user input")
        print("usage: python aggre_stats_mpi.py <Data Path> <Start Date> <End Date> \
												<Polygon boundaries> <Lat & Lon Grid Size > \
												<Sampling number larger than 0> \
												<1/0> <1/0> <1/0> \
												<1/0> <1/0> <1/0> \
												<1/0> <Variable Imput File> <JHist Variable Imput File>")

        # start_date=np.fromstring(sys.argv[2], dtype=np.int, sep='/' )
        # end_date=np.fromstring(sys.argv[2], dtype=np.int, sep='/' )
        # print("Date:",start_date[0],start_date[1],start_date[2])
        # spl_num = np.int(sys.argv[3][1:-1])
        # poly=np.fromstring(sys.argv[1][1:-1], dtype=np.int, sep=',' )
        # grid=np.fromstring(sys.argv[2][1:-1], dtype=np.float, sep=',' )
        # print(spl_num,poly,grid)
        sys.exit()
    else:
        # Define the sampling rate, boundaries of the selected polygon region & the grid size of Lat & Lon
        spl_num = np.int(sys.argv[6][1:-1])
        poly = np.fromstring(sys.argv[4][1:-1], dtype=np.int, sep=',')
        grid = np.fromstring(sys.argv[5][1:-1], dtype=np.float, sep=',')

        # Define the statistics names for HDF5 output
        sts_name = ['Minimum', 'Maximum', 'Mean', 'Pixel_Counts', \
                    'Standard_Deviation', 'Histogram_Counts', 'Jhisto_vs_']

        # Pass system arguments to the function
        sts_switch = np.array(sys.argv[7:14], dtype=np.int)
        sts_switch = np.array((sts_switch == 1))
        varlist = sys.argv[14]

        # Read the variable names from the variable name list
        text_file = np.array(pd.read_csv(varlist, header=0, delim_whitespace=True))  # open(varlist, "r")
        varnames = text_file[:, 0]

        if sts_switch[5] == True:
            intervals_1d = text_file[:, 1]  # This is a string interval arrays
        else:
            intervals_1d = [0]

        if sts_switch[6] == True:
            # Read the joint histogram names from the variable name list
            jvarlist = sys.argv[15]
            text_file = np.array(pd.read_csv(jvarlist, header=0, delim_whitespace=True))  # open(varlist, "r")
            histnames = text_file[:, 1]
            var_idx = text_file[:, 2]  # This is the index of the input variable name which is used for 2D histogram
            intervals_2d = text_file[:, 3]
        else:
            intervals_2d, var_idx = [0], [0]

    # -------------STEP 1: Set up the specific directory --------
    data_path_file = np.array(pd.read_csv(sys.argv[1], header=0, delim_whitespace=True))
    MYD06_dir = data_path_file[0, 0]  # '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD06_L2/'
    MYD06_prefix = data_path_file[0, 1]  # 'MYD06_L2.A'
    MYD03_dir = data_path_file[1, 0]  # '/umbc/xfs1/cybertrn/common/Data/Satellite_Observations/MODIS/MYD03/'
    MYD03_prefix = data_path_file[1, 1]  # 'MYD03.A'
    fileformat = 'hdf'

    # -------------STEP 2: Set up spactial and temporal resolution & variable names----------
    NTA_lats = [poly[0], poly[1]]  # [  0,40] #[-90,90]   #[-30,30]
    NTA_lons = [poly[2], poly[3]]  # [-40,60] #[-180,180] #[-60,60]

    gap_x, gap_y = grid[1], grid[0]  # 0.5,0.625

    if ((NTA_lons[-1] - NTA_lons[0]) % gap_x != 0) | ((NTA_lats[-1] - NTA_lats[0]) % gap_y != 0):
        print("Grid size should be dividable by the dimension of the selected region.")
        print(
            "If you choose the region of latitude  from -40 to 40, then you gird size (Latitude ) should be dividable by 80.")
        print(
            "If you choose the region of longitude from  20 to 35, then you gird size (Longitude) should be dividable by 55.")
        print("Please try again!")
        sys.exit()

    map_lon = np.arange(NTA_lons[0], NTA_lons[1], gap_x)
    map_lat = np.arange(NTA_lats[0], NTA_lats[1], gap_y)
    Lon, Lat = np.meshgrid(map_lon, map_lat)
    grid_lon = np.int((NTA_lons[-1] - NTA_lons[0]) / gap_x)
    grid_lat = np.int((NTA_lats[-1] - NTA_lats[0]) / gap_y)

    # --------------STEP 3: Create arrays for level-3 statistics data-------------------------
    grid_data = {}
    bin_num1 = np.zeros(len(varnames)).astype(np.int)
    bin_num2 = np.zeros(len(varnames)).astype(np.int)
    key_idx = 0
    for key in varnames:
        if sts_switch[0] == True:
            grid_data[key + '_' + sts_name[0]] = np.zeros(grid_lat * grid_lon) + np.inf
        if sts_switch[1] == True:
            grid_data[key + '_' + sts_name[1]] = np.zeros(grid_lat * grid_lon) - np.inf
        if (sts_switch[2] == True) | (sts_switch[3] == True) | (sts_switch[4] == True):
            grid_data[key + '_' + sts_name[2]] = np.zeros(grid_lat * grid_lon)
            grid_data[key + '_' + sts_name[3]] = np.zeros(grid_lat * grid_lon)
            grid_data[key + '_' + sts_name[4]] = np.zeros(grid_lat * grid_lon)
        if sts_switch[5] == True:
            bin_interval1 = np.fromstring(intervals_1d[key_idx], dtype=np.float, sep=',')
            bin_num1[key_idx] = bin_interval1.shape[0] - 1
            grid_data[key + '_' + sts_name[5]] = np.zeros((grid_lat * grid_lon, bin_num1[key_idx]))

            if sts_switch[6] == True:
                bin_interval2 = np.fromstring(intervals_2d[key_idx], dtype=np.float, sep=',')
                bin_num2[key_idx] = bin_interval2.shape[0] - 1
                grid_data[key + '_' + sts_name[6] + histnames[key_idx]] = np.zeros(
                    (grid_lat * grid_lon, bin_num1[key_idx], bin_num2[key_idx]))

        key_idx += 1

    # --------------STEP 4: Read the filename list for different time period-------------------
    fname1, fname2 = [], []

    start_date = np.fromstring(sys.argv[2], dtype=np.int, sep='/')
    end_date = np.fromstring(sys.argv[3], dtype=np.int, sep='/')
    start = date(start_date[0], start_date[1], start_date[2])
    until = date(end_date[0], end_date[1], end_date[2])

    for dt in rrule(DAILY, interval=1, dtstart=start, until=until):
        year = np.int(dt.strftime("%Y"))
        month = np.int(dt.strftime("%m"))
        day = np.int(dt.strftime("%d"))

        data = datetime(year, month, day)
        daynew = data.toordinal()
        yearstart = datetime(year, 1, 1)
        day_yearstart = yearstart.toordinal()
        day_in_year = (daynew - day_yearstart) + 1

        yc = '%04i' % year
        dc = '%03i' % day_in_year

        fname_tmp1 = read_filelist(MYD06_dir, MYD06_prefix, yc, dc, fileformat)
        fname_tmp2 = read_filelist(MYD03_dir, MYD03_prefix, yc, dc, fileformat)
        fname1 = np.append(fname1, fname_tmp1)
        fname2 = np.append(fname2, fname_tmp2)
    print(year, month)

    filenum = np.arange(len(fname1))
    print(len(fname1))

    # --------------STEP 5: Read Attributes of each variables----------------------------------
    unit_list = []
    scale_list = []
    offst_list = []
    longname_list = []
    fillvalue_list = []

    ncfile = Dataset(fname1[0], 'r')

    # Read the User-defined variables from MYD06 product
    tmp_idx = 0
    for key in varnames:
        if key == 'cloud_fraction':
            name_idx = tmp_idx
            continue  # Ignoreing Cloud_Fraction from the input file
        else:
            tmp_data, lonam, unit, fill, scale, offst = readEntry(key, ncfile)
            unit_list = np.append(unit_list, unit)
            scale_list = np.append(scale_list, scale)
            offst_list = np.append(offst_list, offst)
            longname_list = np.append(longname_list, lonam)
            fillvalue_list = np.append(fillvalue_list, fill)
            tmp_idx += 1

    # Add the long name of cloud freaction at the first row
    CM_unit = 'none'
    CM_longname = 'Cloud Fraction from Cloud Mask (cloudy & prob cloudy)'
    CM_fillvalue = -9999
    CM_scale_factor = 0.0001
    CM_add_offset = 0.0
    unit_list = np.insert(unit_list, name_idx, CM_unit)
    scale_list = np.insert(scale_list, name_idx, CM_scale_factor)
    offst_list = np.insert(offst_list, name_idx, CM_add_offset)
    longname_list = np.insert(longname_list, name_idx, CM_longname)
    fillvalue_list = np.insert(fillvalue_list, name_idx, CM_fillvalue)

    ncfile.close()

    # --------------STEP 6: Start Aggregation------------------------------------------------

    # Start counting operation time
    start_time = timeit.default_timer()

    grid_data = run_modis_aggre(fname1, fname2, NTA_lats, NTA_lons, grid_lon, grid_lat, gap_x, gap_y, filenum, \
                                grid_data, sts_switch, varnames, intervals_1d, intervals_2d, var_idx)

    # Compute the mean cloud fraction & Statistics (Include Min & Max & Standard deviation)

    # Reference for statstic parameters
    # sts_name[0]: min
    # sts_name[1]: max
    # sts_name[2]: mean / total_value
    # sts_name[3]: count
    # sts_name[4]: square
    # sts_name[5]: histogram
    # sts_name[6]: joint histogram

    sts_idx = np.array(np.where(sts_switch == True))[0]
    print("Index of User-defined Statistics:", sts_idx)
    key_idx = 0
    for key in varnames:
        for i in sts_idx:
            if i == 0:
                grid_data[key + '_' + sts_name[0]] = grid_data[key + '_' + sts_name[0]].reshape([grid_lat, grid_lon])
            elif i == 1:
                grid_data[key + '_' + sts_name[1]] = grid_data[key + '_' + sts_name[1]].reshape([grid_lat, grid_lon])
            elif i == 2:
                grid_data[key + '_' + sts_name[2]] = (
                        grid_data[key + '_' + sts_name[2]] / grid_data[key + '_' + sts_name[3]])
                grid_data[key + '_' + sts_name[2]] = grid_data[key + '_' + sts_name[2]].reshape([grid_lat, grid_lon])
            elif i == 3:
                grid_data[key + '_' + sts_name[3]] = grid_data[key + '_' + sts_name[3]].reshape([grid_lat, grid_lon])
            elif i == 4:
                grid_data[key + '_' + sts_name[4]] = ((grid_data[key + '_' + sts_name[4]] / grid_data[
                    key + '_' + sts_name[3]].ravel()) - grid_data[key + '_' + sts_name[2]].ravel() ** 2) ** 0.5
                grid_data[key + '_' + sts_name[4]] = grid_data[key + '_' + sts_name[4]].reshape([grid_lat, grid_lon])
            elif i == 5:
                grid_data[key + '_' + sts_name[5]] = grid_data[key + '_' + sts_name[5]].reshape(
                    [grid_lat, grid_lon, bin_num1[key_idx]])
            elif i == 6:
                grid_data[key + '_' + sts_name[6] + histnames[key_idx]] = grid_data[
                    key + '_' + sts_name[6] + histnames[key_idx]].reshape(
                    [grid_lat, grid_lon, bin_num1[key_idx], bin_num2[key_idx]])
        key_idx += 1

    end_time = timeit.default_timer()

    # print('Mean_Fraction:')
    # print( Mean_Fraction  )

    print("Operation Time in {:7.2f} seconds".format(end_time - start_time))

    # --------------STEP 7:  Create HDF5 file to store the result------------------------------
    l3name = 'MYD08_D3' + 'A{:04d}{:02d}'.format(year, month)
    subname = '_baseline_daily_v9_5.h5'
    ff = h5py.File(l3name + subname, 'w')

    PC = ff.create_dataset('lat_bnd', data=map_lat)
    PC.attrs['units'] = 'degrees'
    PC.attrs['long_name'] = 'Latitude_boundaries'

    PC = ff.create_dataset('lon_bnd', data=map_lon)
    PC.attrs['units'] = 'degrees'
    PC.attrs['long_name'] = 'Longitude_boundaries'

    for i in range(sts_idx.shape[0]):
        cnt = 0
        for key in grid_data:

            if key.find("1km") != -1:
                new_name = key.replace("_1km", "")
            else:
                new_name = key

            if (sts_name[sts_idx[i]] in key) == True:
                # print(sts_name[sts_idx[i]],key,grid_data[key].shape)
                # print(longname_list[cnt][:20],new_name)
                addGridEntry(ff, new_name, unit_list[cnt], longname_list[cnt], fillvalue_list[cnt], scale_list[cnt],
                             offst_list[cnt], grid_data[key])
                cnt += 1

    ff.close()

    print(l3name + subname + ' Saved!')
# ---------------------------COMPLETED------------------------------------------------------
