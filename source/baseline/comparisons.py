#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Wed Feb 27 06:17:38 2019
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
"""
import numpy as np
import matplotlib.pyplot as plt
import h5py
from netCDF4 import Dataset
def readData(file):
    f=h5py.File(file,'r')
    CF=f['CF'][:]
    lo=f['lat_bnd'][:]
    la=f['lon_bnd'][:]
    f.close()
    return CF,lo,la
def doPlot(results,benchmark,fig_ttl):
    '''
    results,benchmark:180*360 array
    '''
    X={}
    X['rwDay']=results
    X['bmDay']=benchmark
    fig1,ax1=plt.subplots(3,1,figsize=(6,10))
    fig1.suptitle(fig_ttl)
    cm1=ax1[0].imshow(X['rwDay'],extent=(-180,180,-90,90))
    cm2=ax1[1].imshow(X['bmDay'],extent=(-180,180,-90,90))
    bias=X['rwDay']-X['bmDay']
    cm3=ax1[2].imshow(bias,vmin=-1,vmax=1,extent=(-180,180,-90,90),cmap=plt.cm.seismic)
    ax1[0].set_title('Results')
    ax1[1].set_title('Benchmark')
    ax1[2].set_title('Bias')
    fig1.colorbar(cm1,ax=ax1[0])
    fig1.colorbar(cm2,ax=ax1[1])
    fig1.colorbar(cm3,ax=ax1[2])
    fig1.show()
    return fig1,fig_ttl
if __name__=='__main__':
    sv_results="/umbc/xfs1/jianwu/common/MODIS_Aggregation/output_final2.hdf5"
    zz_results="/umbc/xfs1/jianwu/users/charaj1/CMAC/zz_MODIS_aggregation/Cloud_fraction_test_ZZ.hdf5"
    bmDay_resu="/umbc/xfs1/jianwu/users/charaj1/CMAC/MODIS-Aggregation/MODAgg_3var_oneDay_20080101.hdf5"#single day
    bmTst_resu="/umbc/xfs1/jianwu/users/charaj1/CMAC/MODIS-Aggregation/MODAgg_3var_test_20080101.hdf5"#test (a few files)
    dpDay_resu="/umbc/xfs1/jianwu/common/MODIS_Aggregation/deepak-code/Pandas_one_day_file.hdf"
    rwDay_resu="/umbc/xfs1/jianwu/users/rwalid1/individual/work/Cybtrn-team3/Jan_CF_XR_L2L3_Jun19_n4_run1_1day.nc"
    CF,lo,la={'sv':[],'zz':[],'dpDay':[],'bmDay':[],'rwDay':[]},\
             {'sv':[],'zz':[],'dpDay':[],'bmDay':[],'rwDay':[]},\
             {'sv':[],'zz':[],'dpDay':[],'bmDay':[],'rwDay':[]}
       
    CF['sv'],lo['sv'],la['sv']=readData(sv_results)
    CF['zz'],lo['zz'],la['zz']=readData(zz_results)
    CF['bmDay'],lo['bmDay'],la['bmDay']=readData(bmDay_resu)
#    CF['dp'],lo['dp'],la['dp']=readData(dp_results) #Deepak did not use the latest function.
    f=h5py.File(dpDay_resu,'r')
    CF['dpDay']=f['CM'][:]
    lo['dpDay']=f['Latitude'][:]
    la['dpDay']=f['Longitude'][:]
    f.close()    
    #Redwan has *.nc outputs
    f = Dataset(rwDay_resu, "r")
    CF['rwDay']=f.variables['__xarray_dataarray_variable__'][:]
    f.close()
    fig1,fig1_ttl=doPlot(CF['bmDay'],CF['rwDay'],rwDay_resu.split('/',-1)[-1])

    
#    fig2,ax2=plt.subplots()
#    fig2.suptitle('ZZ_SV_comparison')
#    ax2.plot(CF['bmDay'],CF['rwDay'],'b.')
#    ax2.set_xlabel('ZZ')
#    ax2.set_ylabel('SV')
#    fig2.show()
    
