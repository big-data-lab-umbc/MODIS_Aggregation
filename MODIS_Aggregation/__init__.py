"""Top-level package for MODIS."""

__author__ = """Sukhada Deshpande"""
__email__ = 'sdeshpa1@umbc.edu'
__version__ = '0.1.0'


from .cloud_fraction_aggregate import *
from .baseline_series import *
from .checkaddition import *

# if somebody does "from Sample import *", this is what they will
# be able to access:
__all__ = [
    'aggregateOneFileData'
    ,'displayOutput'
    ,'calculateCloudFraction'
    ,'getInputDirectories'
    ,'read_filelist'
    ,'readEntry'
    ,'read_MODIS'
    ,'cal_stats'
    ,'run_modis_aggre'
    ,'addGridEntry'
    ,'addition'
]
