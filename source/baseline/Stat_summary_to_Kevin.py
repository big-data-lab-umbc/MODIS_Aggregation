#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Mon Sep 30 14:47:59 2019
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
MODIS Aggregation Statistics
Summary to Kevin
**NOT A SCRIPT**
"""

    def Aggregate(self,MOD03_path,MOD06_path,fname_ap=None):
        '''
        - Defines self.l3name: level-3 output file name.
        - Defines self.M (Object) which contains all the computed statistics and variables. 
        fname_ap (string): 'a_string_to_append_to_the_L3_output_file'
        '''
        '''
        For a given MOD03 file path and MOD06 file path
        ***********************************************************************
        This part can be replaced with yours......

        '''
        dt=self.start.split('/')
        yr = [int(dt[2])]
        mn = [int(dt[0])] #np.arange(1,13)  #[1]
        dy = [int(dt[1])] #np.arange(1,32) # [1] #np.arange(1,31)
        # latitude and longtitude boundaries of level-3 grid
        if fname_ap is None:
            self.l3name='MOD08_'+self.l3product+'A{:04d}{:03d}'.format(yr[0],day_of_year(yr[0],mn[0],dy[0]))
        else:
            self.l3name='MOD08_'+self.l3product+'A{:04d}{:03d}'.format(yr[0],day_of_year(yr[0],mn[0],dy[0]))+fname_ap
        lat_bnd = np.arange(-90,91,1)
        lon_bnd = np.arange(-180,180,1)
        nlat = 180
        nlon = 360
    
        '''
        Initializations (Both variables and functions for multiple statistics)
        ***************************************************************************
        Memory allocation and computations are only be done for the requested variables
        (any combination of variables ex. CTP,CTT,COT,CER,etc.) and statistics
        (any combination of min,max,mean,stdd).
        
        o We don't have to deal with different combinations. So a single function to compute all the stats.
        '''    
        M=Memory() # An empty object to store variables
        M.TOT_pix      = np.zeros(nlat*nlon)#To compute CF 
        M.CLD_pix      = np.zeros(nlat*nlon)#Will be needed for all others including CF
        M.CF = {'CF':np.zeros(nlat*nlon),'min':np.zeros(nlat*nlon),'max':np.zeros(nlat*nlon)}
        
        variable_names={'cloud_optical_thickness_1km','cloud_top_temperature_1km'} #list of level 2 varaibles
        M.stt['min':{},'max':{},'mean':{},'stdd':{}]#dictionary will containes 'variable_names'. 'mean' will be needed for stdd calculations.
        #mean and stdd (Initialization)
        M.XXX_pix={}
        for key in variable_names:
            M.XXX_pix[key]=np.zeros(nlat*nlon)
        M.XXX_pixSq=M.XXX_pix # For Stdd
        #Min and Max (Initialization) 
        M.mnx={};M.stt={}
        M.mnx['min']={}
        for key in variable_names:
            M.mnx['min'][key]=np.zeros(nlat*nlon)+np.inf
        M.mnx['max']={}
        for key in self.variables:
            M.mnx['max'][key]=np.zeros(nlat*nlon)-np.inf
        '''
        ***********************************************************************
        ***********************************************************************
        '''
        # Min and Max computations (function definitions)

            #Both min and max are needed
            def minmax(val,j,M):
                mn,mx=val.min(),val.max()
                if mn<M.mnx['min'][key][j]:
                    M.mnx['min'][key][j]=mn
                if mx>M.mnx['max'][key][j]:
                    M.mnx['max'][key][j]=mx

        '''
        ***********************************************************************
        '''        
            #Mean and std. computations. minmax() is called inside this function
            M.stt['mean'],M.stt['stdd']={},{}# Mean is needed to calculate Std
            def MeanStd(data,j,latlon_index,M):
                #Both mean and stdd
                for key in data:
                    #print(key)
                    val=data[key][np.where(latlon_index == j)]
                    M.XXX_pix[key][j]=M.XXX_pix[key][j]+np.sum(val)
                    M.XXX_pixSq[key][j]=M.XXX_pixSq[key][j]+np.sum(val**2)   
                    minmax(val,j,M)

        '''
        ***********************************************************************
        ***********************************************************************
        '''
        '''
        Looping through files
        ======================
        '''
        #-----------------------------------------------
        tot_F=0 #Total number of file couples read
        self.start=time.time()
        MD=MODIS_Level2(self.variables)
        for y,m,d in  itertools.product(yr,mn, dy):
            #-------------find the MODIS prodcts--------------#
            #date = datetime.datetime(y,m,d)
            JD01, JD02 = gcal2jd(y,1,1)
            JD1, JD2 = gcal2jd(y,m,d)
            JD = np.int((JD2+JD1)-(JD01+JD02) + 1)
            granule_time = datetime.datetime(y,m,d,0,0)
            while granule_time <= datetime.datetime(y,m,d,23,55):  # 23,55
    #            print('granule time:',granule_time)
                MOD03_fp = 'MYD03.A{:04d}{:03d}.{:02d}{:02d}.006.?????????????.hdf'.format(y,JD,granule_time.hour,granule_time.minute)
                MOD06_fp = 'MYD06_L2.A{:04d}{:03d}.{:02d}{:02d}.006.?????????????.hdf'.format(y,JD,granule_time.hour,granule_time.minute)
                MOD03_fn, MOD06_fn =[],[]
                for MOD06_flist in  os.listdir(MOD06_path):
                    if fnmatch.fnmatch(MOD06_flist, MOD06_fp):
                        MOD06_fn = MOD06_flist
                for MOD03_flist in  os.listdir(MOD03_path):
                    if fnmatch.fnmatch(MOD03_flist, MOD03_fp):
                        MOD03_fn = MOD03_flist
                if MOD03_fn and MOD06_fn: # if both MOD06 and MOD03 products are in the directory
                    tot_F+=1
                    Lat,Lon,data = MD.readHDF(MOD03_path+MOD03_fn,MOD06_path+MOD06_fn)
                    CM=data['CM'].ravel()
                    del data['CM'] #CM is not a final product
                    for key in data:
                        data[key]=data[key].ravel()    
                    Lat=Lat.ravel()
                    Lon=Lon.ravel()
                    lat_index = value_locate(lat_bnd,Lat)
                    lon_index = value_locate(lon_bnd,Lon)
                    latlon_index = lat_index*nlon + lon_index
                    latlon_index_unique = np.unique(latlon_index)
                    
                    for i in np.arange(latlon_index_unique.size):
                        j=latlon_index_unique[i]
                        M.TOT_pix[j] = M.TOT_pix[j]+np.sum(CM[np.where(latlon_index == j)]>=0)
                        M.CLD_pix[j] = M.CLD_pix[j]+np.sum(CM[np.where(latlon_index == j)]<=1)
                        #To calculate other variables and statistics---------------------------
                        '''
                        *********************************************************************
                        Here we call MeanStd()
                        '''
                        MeanStd(data,j,latlon_index,M)#Here we call MeanStd function

                    

                granule_time += datetime.timedelta(minutes=5)
        '''
        ***************************************************************
        Computing all the statistics
        '''
        #Cloud fractions
        M.total_cloud_fraction  =  division(M.CLD_pix,M.TOT_pix).reshape([nlat,nlon])
        M.pixel_count = M.CLD_pix.reshape([nlat,nlon])
        #The other statistics
        for key in data:
            for st in M.stt:
                if st == 'stdd':
                    M.stt['mean'][key]=division(M.XXX_pix[key],M.CLD_pix).reshape([nlat,nlon])
                    #stdd=np.sqrt(<Xi^2>-<X>^2)
                    M.stt[st][key]=np.sqrt(division(M.XXX_pixSq[key],M.CLD_pix).reshape([nlat,nlon])-M.stt['mean'][key]**2)
                elif st == 'mean':
                    M.stt[st][key]=division(M.XXX_pix[key],M.CLD_pix).reshape([nlat,nlon])
                if st == 'min' or st == 'max':
                    M.stt[st][key]=M.mnx[st][key].reshape([nlat,nlon])
        '''
        M contains all the stuff
        '''
        
        
