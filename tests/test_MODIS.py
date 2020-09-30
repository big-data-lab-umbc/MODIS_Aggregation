import unittest
from MODIS_Aggregation import getInputDirectories
from MODIS_Aggregation import aggregateOneFileData
from MODIS_Aggregation import displayOutput

class getFilePathTest(unittest.TestCase):

    def test_valid_file_path(self):
        x, y = getInputDirectories()
        self.assertIsNotNone(x)
        self.assertIsNotNone(y)

class TruthTest(unittest.TestCase):

    def test_assert_true(self):
            self.assertTrue(True)

    def test_assert_false(self):
            self.assertFalse(False)

    # def test_invalid_file_path(self):
    #     x, y = getInputDirectories()
    #     self.assertIsNone(x)
    #     self.assertIsNone(y)

    # def test_valid_aggregate_data(self):
    #     a,b = aggregateOneFileData()
    #     self.assertIsNotNone(a)
    #     self.assertIsNotNone(b)
    #
    # def test_display_output(self):
    #     p = displayOutput(self)
    #     self.assertIsNotNone(displayOutput(cf))

if __name__ == '__main__':
    unittest.main()

# #!/usr/bin/env python
#
# """Tests for `MODIS_Aggregation` package."""
#
#
# import unittest
# from click.testing import CliRunner
#
# from MODIS_Aggregation import *
# from MODIS_Aggregation import cli
#
#
# class TestModis(unittest.TestCase):
#     """Tests for `MODIS_Aggregation` package."""
#
#     def setUp(self):
#         """Set up test fixtures, if any."""
#
#     def tearDown(self):
#         """Tear down test fixtures, if any."""
#
#     def test_000_something(self):
#         """Test something."""
#
#     def test_command_line_interface(self):
#         """Test the CLI."""
#         runner = CliRunner()
#         result = runner.invoke(cli.main)
#         assert result.exit_code == 0
#         assert 'MODIS.cli.main' in result.output
#         help_result = runner.invoke(cli.main, ['--help'])
#         assert help_result.exit_code == 0
#         assert '--help  Show this message and exit.' in help_result.output


