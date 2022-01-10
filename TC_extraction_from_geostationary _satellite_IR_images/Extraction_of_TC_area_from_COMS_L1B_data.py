#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 30 10:46:22 2021

@author: yhbaek
"""


import sys
import os
import time as TT
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt

# %%

IBTrACS_path='/home/yhbaek/STUDY/DATA/Best-Track/IBTrACS/IBTrACS_JTWC_20210822.npy'
IBTrACS=np.load(IBTrACS_path)

find_exist_data=np.where(np.isnan(IBTrACS[:, :, 0])==0)


New_IBTrACS=np.zeros( (len(find_exist_data[0]), IBTrACS.shape[2]) )

for L in range(len(find_exist_data[0])):
    New_IBTrACS[L, :]=IBTrACS[find_exist_data[0][L], find_exist_data[1][L], :]
del(L)    
del(find_exist_data, IBTrACS, IBTrACS_path)    

New_IBTrACS=New_IBTrACS[np.where( (np.isnan(New_IBTrACS[:, 6])==0) & (np.isnan(New_IBTrACS[:, 7])==0))[0], :]

New_IBTrACS=New_IBTrACS[np.where( (New_IBTrACS[:, 2]>=2011) & (New_IBTrACS[:, 2]<=2020) )[0], :]

find_COMS_time=np.where( (New_IBTrACS[:, 5]==2) |  (New_IBTrACS[:, 5]==5) | (New_IBTrACS[:, 5]==8) | 
                        (New_IBTrACS[:, 5]==11) |  (New_IBTrACS[:, 5]==14) | (New_IBTrACS[:, 5]==17) | 
                        (New_IBTrACS[:, 5]==20) |  (New_IBTrACS[:, 5]==23) )[0]

New_IBTrACS=New_IBTrACS[find_COMS_time, :]

del(find_COMS_time)
# %%
DIR='/home/data/COMS/L1B/'
Pixels=300 # 300pixels*4km = 1200km * 2 = 2400km

COMS_Lat=np.load(DIR + 'COMS_lon_lat_infomation/COMS_Lat.npy')
COMS_Lon=np.load(DIR + 'COMS_lon_lat_infomation/COMS_Lon.npy')


missing_data=np.zeros((New_IBTrACS.shape[0], 5))
ct=0
#L=4378
for L in range(New_IBTrACS.shape[0]):
#for L in np.arange(4378, 14170):    
    Year=str(int(New_IBTrACS[L, 2]))
    Month=str(int(New_IBTrACS[L, 3])).zfill(2)
    Day=str(int(New_IBTrACS[L, 4])).zfill(2)
    Hour=str(int(New_IBTrACS[L, 5])).zfill(2)
    
    Date=Year + Month + Day

    Out_DIR=DIR + 'TC_Area/' + Year + '/' + Month

    os.makedirs(Out_DIR, exist_ok=True)        

    COMS_path=DIR + 'FD/' + Year + Month + '/' + Day  + '/coms_mi_le1b_all_fd000ge_' + Date + Hour + '15.he5'

    
    if os.path.exists(COMS_path)==True:
        
        COMS_file=Dataset(COMS_path, 'r')
        COMS_group=COMS_file.groups['HDFEOS']['GRIDS']['IR Image Pixel Value']['Data Fields']

        # %% Cut select TC Area
        TC_Lon=New_IBTrACS[L, 6]
        TC_Lat=New_IBTrACS[L, 7]
        
        diff=np.sqrt( ((COMS_Lon-TC_Lon)**2) + ((COMS_Lat-TC_Lat)**2) )
        
        x=np.where(diff==np.nanmin(diff))[0][0]
        y=np.where(diff==np.nanmin(diff))[1][0]
        
        x1=x-Pixels
        x2=x+Pixels+1
        
        y1=y-Pixels
        y2=y+Pixels+1
        

        if (x1>=0) & (y1>=0) & (x2<=COMS_Lon.shape[0]) & (y2<=COMS_Lat.shape[0]):
            cut_lat=COMS_Lat[x1:x2, y1:y2]
            cut_lon=COMS_Lon[x1:x2, y1:y2]
        
            IR1=np.array(COMS_group.variables['IR1 band Image Pixel Values'][x1:x2, y1:y2], dtype=float)
            IR2=np.array(COMS_group.variables['IR2 band Image Pixel Values'][x1:x2, y1:y2], dtype=float)
            SW=np.array(COMS_group.variables['SWIR band Image Pixel Values'][x1:x2, y1:y2], dtype=float)
            WV=np.array(COMS_group.variables['WV band Image Pixel Values'][x1:x2, y1:y2], dtype=float)
            
            #IR1[np.where(IR1==32768)]=np.NaN
            #IR2[np.where(IR2==32768)]=np.NaN
            #SW[np.where(SW==32768)]=np.NaN
            #WV[np.where(WV==32768)]=np.NaN
                    
            del(COMS_file, COMS_group)
            
            # %% 
            COMS_table=pd.read_csv(DIR + "Conversion_Table/COMS_MI_Conversion_Table_new_20111201.csv")
            
            Digital_Count=np.array(COMS_table['Digital_Count'])
            
            IR1_table=np.array(COMS_table['IR1_TBB'])
            IR2_table=np.array(COMS_table['IR2_TBB'])
            SW_table=np.array(COMS_table['SWIR_TBB'])
            WV_table=np.array(COMS_table['WV_TBB'])
            
                
            for D in range(IR1_table.shape[0]):
                IR1[np.where(IR1==D)]=IR1_table[D]
                IR2[np.where(IR2==D)]=IR2_table[D]
                SW[np.where(SW==D)]=SW_table[D]
                WV[np.where(WV==D)]=WV_table[D]
            del(D)
            
            """
            Nor_IR1=IR1/(np.where(IR1_table==0)[0][0]-1)            # 950
            Nor_IR2=IR2/(np.where(IR2_table==0)[0][0]-1)            # 951
            Nor_SW=SW/(np.where(SW_table==0)[0][0]-1)               # 950
            Nor_WV=WV/(np.where(WV_table==0)[0][0]-1)               # 949
            """
        
            del(COMS_table, Digital_Count)
            del(IR1_table, IR2_table, SW_table, WV_table)
        
            # %% Make output
            TC_info=np.zeros((25, 2), dtype='<U13')
                    
            TC_info[:, 0]=np.array(['Storm ID', 'Date', 'Hour', 'Longitude ', 'Latitude', 
                 'MSW', 'DIST2LAND', 'Pres', 'POCI', 'EYE', 'RMW', 'ROCI',
                 'R34_Q1', 'R34_Q2', 'R34_Q3', 'R34_Q4',
                 'R50_Q1', 'R50_Q2', 'R50_Q3', 'R50_Q4',
                 'R64_Q1', 'R64_Q2', 'R64_Q3', 'R64_Q4',
                 'Heading_angle'])  
    
            TC_info[0, 1]=int(New_IBTrACS[L, 0])
            TC_info[1, 1]=Date
            TC_info[2, 1]=Hour
            TC_info[3:, 1]=New_IBTrACS[L, 6:]
    
            
            NCFILE=Dataset(Out_DIR + '/COMS_MI_LE1B_TC_AREA_040GE_rot_000_' + Date + Hour + '00_' + str(int(New_IBTrACS[L, 0])) + '.nc',
                           'w', format='NETCDF4')
            
            info_i=NCFILE.createDimension('info_i', TC_info.shape[0])
            info_j=NCFILE.createDimension('info_j', TC_info.shape[1])
            lat=NCFILE.createDimension('Latitude', cut_lat.shape[0])
            lon=NCFILE.createDimension('Longitude', cut_lon.shape[0])
            
            lats=NCFILE.createVariable('lat', 'f8', ('Longitude', 'Latitude'))
            lons=NCFILE.createVariable('lon', 'f8', ('Longitude', 'Latitude'))
    
            ir1=NCFILE.createVariable('IR1', 'f8', ('Longitude', 'Latitude'))
            ir2=NCFILE.createVariable('IR2', 'f8', ('Longitude', 'Latitude'))
            sw=NCFILE.createVariable('SW', 'f8', ('Longitude', 'Latitude'))
            wv=NCFILE.createVariable('WV', 'f8', ('Longitude', 'Latitude'))
    
            #nor_ir1=NCFILE.createVariable('Nor_IR1', 'f8', ('Longitude', 'Latitude'))
            #nor_ir2=NCFILE.createVariable('Nor_IR2', 'f8', ('Longitude', 'Latitude'))
            #nor_sw=NCFILE.createVariable('Nor_SW', 'f8', ('Longitude', 'Latitude'))
            #nor_wv=NCFILE.createVariable('Nor_WV', 'f8', ('Longitude', 'Latitude'))
            
            tc_info=NCFILE.createVariable('TC_info', 'S8', ('info_i', 'info_j'))
    
            NCFILE.description = 'Extract only TC area from COMS MI LE1B 4km data'
            NCFILE.history = 'Created on ' + TT.ctime(TT.time())
            NCFILE.author = 'Y-H Baek, Typhoon research center, Jeju national university'
            lats.units = 'degrees'
            lons.units = 'degrees'
            ir1.units = 'IR1 Brightness temperature (TBB), unit:K'
            ir2.units = 'IR2 Brightness temperature (TBB), unit:K'
            sw.units = 'SW Brightness temperature (TBB), unit:K'
            wv.units = 'WV Brightness temperature (TBB), unit:K'
    
            #nor_ir1.units = 'Normalize by dividing by the largest of the IR1 correction values(950)'
            #nor_ir2.units = 'Normalize by dividing by the largest of the IR2 correction values(951)'
            #nor_sw.units = 'Normalize by dividing by the largest of the SW correction values(950)'
            #nor_wv.units = 'Normalize by dividing by the largest of the WV correction values(949)'
            
            tc_info.units = 'MWS : kt, DIST2LAND : km, Pres : hPa, POCI : hPa, EYE : n mi, RMW : n mi, ROCI : n mi, R34 : n mi, R50 : n mi, R64 : n mi, Heading angle : degrees'
    
            lats[:]=cut_lat
            lons[:]=cut_lon
            
            ir1[:]=IR1
            ir2[:]=IR2
            sw[:]=SW
            wv[:]=WV
            
            #nor_ir1[:]=Nor_IR1
            #nor_ir2[:]=Nor_IR2
            #nor_sw[:]=Nor_SW
            #nor_wv[:]=Nor_WV
            
            tc_info[:]=TC_info
            
            NCFILE.close()
            del(TC_info, info_i, info_j, lat, lon)
            del(lats, lons, ir1, ir2, sw, wv)
            del(tc_info)
            #del(nor_ir1, nor_ir2, nor_sw, nor_wv)
        
            del(cut_lat, cut_lon)
            del(IR1, IR2, SW, WV)
            #del(Nor_IR1, Nor_IR2, Nor_SW, Nor_WV)
            del(diff, x, y, x1, y1, x2, y2, TC_Lon, TC_Lat)
        else:
            missing_data[ct, 0]=Year
            missing_data[ct, 1]=Month
            missing_data[ct, 2]=Day
            missing_data[ct, 3]=Hour
            ct=ct+1

    # %%
    else:
        missing_data[ct, 0]=Year
        missing_data[ct, 1]=Month
        missing_data[ct, 2]=Day
        missing_data[ct, 3]=Hour
        ct=ct+1

    del(COMS_path)

    print(L, '/', str(New_IBTrACS.shape[0]), Date, '-', Hour)
    
    del(Year, Month, Day, Hour, Date)

del(L)
del(COMS_Lon, COMS_Lat, Pixels, New_IBTrACS, ct)
 
missing_data=missing_data[np.where(missing_data[:, 0]!=0)[0], :]
np.save(DIR + 'missing_data_COMS_JTWC.npy', missing_data)
#del(missing_data, Out_DIR)









# %% make GK2A Longitude and Latitude
"""
Resolution=2.0
i=np.arange(0, IR105_file.getncattr('number_of_columns'), dtype='f')
j=np.arange(0, IR105_file.getncattr('number_of_lines'), dtype='f')

i,j=np.meshgrid(i, j)
(GK2A_Lat, GK2A_Lon)=latlon_from_lincol_geos(Resolution, j, i)
del(i, j)

np.save(DIR + 'latlon/GK2A_Lon_2km.npy', GK2A_Lon)
np.save(DIR + 'latlon/GK2A_Lat_2km.npy', GK2A_Lat)
"""




