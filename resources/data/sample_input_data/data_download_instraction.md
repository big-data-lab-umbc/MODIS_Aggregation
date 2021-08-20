# Description
This document is to guide users to download the required MODIS Level-2 data for the example run 

## Requirement of Input Level-2 Data 
Our aggreagtion code requires the input Level-2 data to cover the time range of at least one day with additional 3 hours. The one-day data is minimum for aggregation of daily Level-3 data. The additional 3 hours after the selected day is for the adjustment of the new Definition of Day (Reference: Page 5 to 7 in https://atmosphere-imager.gsfc.nasa.gov/sites/default/files/ModAtmo/documents/L3_ATBD_C6_C61_2020_08_06.pdf).

Note that the wget command for MODIS data download does not allow to download part of data within a day. Therefore, the download process covers two days Level-2 data including MYD06_L2 (Cloud Product) and MYD03 (Geolocation Product), which would take 47G disk space in total. 

## Level-2 Data Download Instrcution (wget)

### Step 1: Create an account or sign in on NASA Earthdata (https://urs.earthdata.nasa.gov/) 

### Step 2: Request a token for using the wget command-line utility in your preferred terminal
In order to properly authenticate your transfer and download, please obtain an app key (https://ladsweb.modaps.eosdis.nasa.gov/tools-and-services/data-download-scripts/#requesting) according to these instructions. (https://ladsweb.modaps.eosdis.nasa.gov/tools-and-services/data-download-scripts/#appkeys) 

### Step 3: Modify the token in the wget command-line of the download_modis.sh 
```
>> vi download_modis.sh
Replace <put your valid download token here> by your token from your account

>> ./download_modis.sh
```
### Step 4: check the [MYD06_L2] & [MYD03] folder foe the downloaded data.

