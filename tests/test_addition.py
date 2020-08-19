import unittest
from MODIS_Aggregation import checkaddition


class CheckAdd(unittest.TestCase):

    def test_addition(self):
        result = checkaddition.addition(3, 2)
        expected = 5
        self.assertEqual(expected, result, msg=None)
        print(result)


if __name__ == '__main__':
    unittest.main()

