#!/bin/bash

echo 'READY TO USE WGET DOWNLOAD FILES!!'
echo 'DOWNLOADING TWO DAYS DATA'

for j in {2008..2008}
do

	for i in {1..2}
	do
	
		wget -c --retry-connrefused --tries=0 --timeout=20 --waitretry=20 -e robots=off -m -np -nd -A hdf -R .html,.tmp -nH --cut-dirs=3 "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/MYD06_L2/${j}/00${i}/" --header "Authorization: Bearer <put your valid download token here>" -P MYD06_L2/
	
	done 
	
	for i in {1..2}
	do
	
		wget -c --retry-connrefused --tries=0 --timeout=20 --waitretry=20 -e robots=off -m -np -nd -A hdf -R .html,.tmp -nH --cut-dirs=3 "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/MYD03/${j}/00${i}/" --header "Authorization: Bearer <put your valid download token here>" -P MYD03/
	
	done
	
done


