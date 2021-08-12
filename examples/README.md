# Description
This folder shows examples on how to use our software.

# Local Execution
Folder local_execution shows our [Dask](https://dask.org/) based sample code to call our MODIS aggregation code to run serially (no parallelization) on a single machine. MODIS_Aggregation_Local.py will call MODIS_Aggregation functions to conduct aggregation.

# Dask based Distributed Execution
Folder dask_based_distributed_execution shows our [Dask](https://dask.org/) based sample code to call our MODIS aggregation code to run on a distributed cluster. MODIS_Aggregation_DASK.slurm will need to be submitted to the cluster's scheduler. Within the slurm file, MODIS_Aggregation_DASK.py is called to call MODIS_Aggregation functions in parallel.

To run the code, additional libraries in the requirements_dask.txt should be installed first.

# MPI based Distributed Execution
Folder mpi_based_distributed_execution shows our [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) based sample code to call our MODIS aggregation code to run on a distributed cluster. MODIS_Aggregation_MPI.slurm will need to be submitted to the cluster's scheduler. Within the slurm file, MODIS_Aggregation_MPI.py is called to call MODIS_Aggregation functions in parallel.

To run the code, additional libraries in the requirements_mpi.txt should be installed first. 

# Auxiliary Files
The csv files on the top folder are example input files the above three execution approaches will use as inputs.

# Result Visualization and Comparison
Folder result_comparison shows our sample Jupyter Notebooks to visualize our results and compare them with MODIS Level 3 products.

# Deployment as Service
Because the deployment involves additional packages, we please check [this document](https://github.com/big-data-lab-umbc/stratus/blob/master/Stratus-MODIS-Aggregation-Deployment.md) at another repository on how to utilize [stratus](https://github.com/nasa-nccs-cds/stratus) to deploy our MODIS aggregation codes as services.
