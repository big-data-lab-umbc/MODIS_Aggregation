# Description
This repository does NASA MODIS data aggregation from level 2 to level 3 in a flexible and parallel approach.
Please check documents in examples folder to see examples on how to use our software.

### Installation
#### Conda environment setup
```
conda create -n MODIS_Aggregation -c conda-forge python=3.7 libnetcdf netCDF4 h5py

>> git clone https://github.com/big-data-lab-umbc/MODIS_Aggregation.git
>> cd MODIS_Aggregation
>> python setup.py install
```

The code is tested with Python 3.7

# Team members
- Faculty: [Dr. Jianwu Wang](https://userpages.umbc.edu/~jianwu/), Department of Information Systems, UMBC
- Faculty: [Dr. Zhibo Zhang](https://physics.umbc.edu/people/faculty/zhang/), Department of Physics, UMBC
- PhD student: Jianyu Zheng, Department of Physics, UMBC
- PhD student: Chamara Rajapakshe, Department of Physics, UMBC
- PhD student: Pei Guo, Department of Information Systems, UMBC
- PhD student: Xin Huang, Department of Information Systems, UMBC
- MS student: Supriya Sangondimath, Department of Information Systems, UMBC
- MS student: [Savio Kay](https://saviokay.com), Department of Information Systems, UMBC
- MS student: Deepak Prakash, Department of Information Systems, UMBC
- MS student: Lakshmi Priyanka Kandoor, Department of Information Systems, UMBC

# Publications
- Jianwu Wang, Xin Huang, Jianyu Zheng, Chamara Rajapakshe, Savio Kay, Lakshmi Kandoor, Thomas Maxwell, Zhibo Zhang. Scalable Aggregation Service for Satellite Remote Sensing Data. In Proceedings of the 20th International Conference on Algorithms and Architectures for Parallel Processing (ICA3PP 2020), pages 184-199, Springer, 2020.

# Acknowledgement
The project is mainly funded by NASA CMAC program
