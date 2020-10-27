import unittest
from MODIS_Aggregation import *
import glob
import time
import xarray as xr
import netCDF4


class GetOutput(unittest.TestCase):


    def test_check_value(self):
        M03_dir = "../resources/data/sample_input_data/MYD03/"
        M06_dir = "../resources/data/sample_input_data/MYD06_L2/"
        print(M06_dir)
        print(M03_dir)
        M03_files = sorted(glob.glob(M03_dir + "MYD03.A2008*"))
        M06_files = sorted(glob.glob(M06_dir + "MYD06_L2.A2008*"))
        print(M03_files)
        print(M06_files)
        t0 = time.time()
        # calculate cloud fraction
        cf = calculateCloudFraction(M03_files, M06_files)
                # calculate execution time
        t1 = time.time()
        total = t1 - t0
        print("total execution time (Seconds):" + str(total))
                # display the output
        result = cf[127, 196]
        print("result:", result)

        expected = 0.12883436
        self.assertAlmostEqual(expected, result, places=8, msg=None)
        print(result)


if __name__ == '__main__':
    unittest.main()