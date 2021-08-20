# Description
This folder shows examples on how to use our software.

# Input Level-2 Data Instruction
This instrcution is to guide users to download the required MODIS Level-2 data for the example run 

## Requirement of Input Level-2 Data 
Our aggreagtion code requires the input Level-2 data to cover the time range of at least one day with additional 3 hours. The one-day data is minimum for aggregation of daily Level-3 data. The additional 3 hours after the selected day is for the adjustment of the new Definition of Day (Reference: Page 5 to 7 in https://atmosphere-imager.gsfc.nasa.gov/sites/default/files/ModAtmo/documents/L3_ATBD_C6_C61_2020_08_06.pdf).

Note that the wget command for MODIS data download does not allow to download part of data within a day. Therefore, the download process covers two days Level-2 data including MYD06_L2 (Cloud Product) and MYD03 (Geolocation Product), which will take 47G disk space in total. 

## Level-2 Data Download Instrcution (wget)

### Step 1: Create an account or sign in on NASA Earthdata (https://urs.earthdata.nasa.gov/) 

### Step 2: Request a token for using the wget command-line utility in your preferred terminal
In order to properly authenticate your transfer and download, please obtain an app key (https://ladsweb.modaps.eosdis.nasa.gov/tools-and-services/data-download-scripts/#requesting) according to these instructions. (https://ladsweb.modaps.eosdis.nasa.gov/tools-and-services/data-download-scripts/#appkeys) 

### Step 3: Modify the token in the wget command-line of the download_modis.sh 
```
>> vi download_modis.sh
Replace <put your valid download token here> by your newly requested token

>> ./download_modis.sh
```
### Step 4: check the MYD06_L2 & MYD03 folder in resources/data/sample_input_data/ for the downloaded data.

## Replacement for other MODIS Level-2 dataset
Users can explore other MODIS Level-2 dataset through https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/. The wget command should be changed accordingly if users want to download MODIS Level-2 dataset other than MYD06_L2 and MYD03.

# Local Execution
Folder local_execution shows our sample code to call our MODIS aggregation code to run serially (no parallelization) on a single machine. MODIS_Aggregation_Local.py will call MODIS_Aggregation functions to conduct aggregation.

# Dask based Distributed Execution
Folder dask_based_distributed_execution shows our [Dask](https://dask.org/) based sample code to call our MODIS aggregation code to run on a distributed cluster. MODIS_Aggregation_DASK.slurm will need to be submitted to the cluster's scheduler to allocate one node. Within the slurm file, MODIS_Aggregation_DASK.py first requests additional compute nodes via its SLURMCluster() and scale() functions, then calls MODIS_Aggregation functions in parallel on the nodes.

To run the code, additional libraries in the requirements_dask.txt should be installed first.

# MPI based Distributed Execution
Folder mpi_based_distributed_execution shows our [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) based sample code to call our MODIS aggregation code to run on a distributed cluster. MODIS_Aggregation_MPI.slurm will need to be submitted to the cluster's scheduler to allocate all nodes required by the execution. Within the slurm file, MODIS_Aggregation_MPI.py calls MODIS_Aggregation functions in parallel on the nodes allocated.

To run the code, additional libraries in the requirements_mpi.txt should be installed first.

# Auxiliary Files
The csv files on the top folder are example input files the above three execution approaches will use as inputs.

# Result Visualization and Comparison
Folder result_comparison shows our sample Jupyter Notebooks to visualize our results and compare them with MODIS Level 3 products.

# Deployment as Service
Because the deployment involves additional packages, we please check [this document](https://github.com/big-data-lab-umbc/stratus/blob/master/Stratus-MODIS-Aggregation-Deployment.md) at another repository on how to utilize [stratus](https://github.com/nasa-nccs-cds/stratus) to deploy our MODIS aggregation codes as services.
