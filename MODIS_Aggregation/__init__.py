"""Top-level package for MODIS."""

__author__ = """X. Huang"""
__email__ = 'xinh1@umbc.edu'
__version__ = '0.1.0'


#from .cfa import *
# from .baseline_series import *
#from .checkaddition import *

from .aggregate_functions import *

# if somebody does "from Sample import *", this is what they will
# be able to access:
__all__ = [
    'read_filelist'
    ,'read_user_inputs'
    ,'readEntry'
    ,'read_MODIS'
    ,'cal_stats'
    ,'run_modis_aggre'
    ,'addGridEntry'
]
