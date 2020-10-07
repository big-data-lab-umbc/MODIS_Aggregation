import unittest
from MODIS_Aggregation import new_func


class GetOutput(unittest.TestCase):

    def test_check_value(self):
        def cfa_check():
            M03_dir, M06_dir = getInputDirectories()
            print(M06_dir)
            print(M03_dir)
            M03_files = sorted(glob.glob(M03_dir + "MYD03.A2008*"))
            M06_files = sorted(glob.glob(M06_dir + "MYD06_L2.A2008*"))
            t0 = time.time()
            # calculate cloud fraction
            cf = calculateCloudFraction(M03_files, M06_files)
            # calculate execution time
            t1 = time.time()
            total = t1 - t0
            print("total execution time (Seconds):" + str(total))
            # display the output
            displayOutput(xr.DataArray(cf))
            return (xr.open_dataset("../tests/monthlyCloudFraction-file-level-for-loop.nc")[
                        '__xarray_dataarray_variable__'][137, 201].values)

        result = new_func.cfa_check()
        expected = 0.0
        self.assertEqual(expected, result, msg=None)
        print(result)


if __name__ == '__main__':
    unittest.main()