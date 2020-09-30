import unittest
from MODIS_Aggregation import new_func


class GetOutput(unittest.TestCase):

    def test_check_value(self):
        result = new_func.cfa_check()
        expected = 0.0
        self.assertEqual(expected, result, msg=None)
        print(result)


if __name__ == '__main__':
    unittest.main()
