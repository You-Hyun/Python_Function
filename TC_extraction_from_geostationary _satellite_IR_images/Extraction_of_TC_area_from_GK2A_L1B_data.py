#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# %%
"""
Created on Thu Aug  5 15:05:52 2021

@author: yhbaek
"""

import sys
import os
import time as TT
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt

#sys.path.append('/data/GK2A/L1B/util/')
#from GK2A_util import cut_with_latlon_geos, latlon_from_lincol_geos, lincol_from_latlon_geos

# %%

IBTrACS_path='/home/yhbaek/STUDY/DATA/Best-Track/IBTrACS/IBTrACS_20211231.npy'

IBTrACS=np.load(IBTrACS_path)

find_exist_data=np.where(np.isnan(IBTrACS[:, :, 0])==0)


New_IBTrACS=np.zeros( (len(find_exist_data[0]), IBTrACS.shape[2]) )

for L in range(len(find_exist_data[0])):
    New_IBTrACS[L, :]=IBTrACS[find_exist_data[0][L], find_exist_data[1][L], :]
del(L)    
del(find_exist_data, IBTrACS, IBTrACS_path)    

New_IBTrACS=New_IBTrACS[np.where( (np.isnan(New_IBTrACS[:, 6])==0) & (np.isnan(New_IBTrACS[:, 7])==0))[0], :]

# %%
# GK2A 2019.07.01 ~ (IOT : 2019.07.01~2020.03.31)

find_COMS_time=np.where( (New_IBTrACS[:, 5]==2) |  (New_IBTrACS[:, 5]==5) | (New_IBTrACS[:, 5]==8) | 
                        (New_IBTrACS[:, 5]==11) |  (New_IBTrACS[:, 5]==14) | (New_IBTrACS[:, 5]==17) | 
                        (New_IBTrACS[:, 5]==20) |  (New_IBTrACS[:, 5]==23) )[0]
New_IBTrACS=New_IBTrACS[find_COMS_time, :]
del(find_COMS_time)


# GK2A 2019.07.01 ~ (IOT : 2019.07.01~2020.03.31)
New_IBTrACS=New_IBTrACS[np.min(np.where( (New_IBTrACS[:, 2]>=2019) & (New_IBTrACS[:, 3]==7))[0]):, :]        


New_IBTrACS=New_IBTrACS[506:1608, :] # IOT IBTrACS



# %%
DIR='/data/GK2A/'
Pixels=600 # 600pixels * 2km(resolution) = 1200km(radius from TC center) * 2 = 2400km

GK2A_Lon=np.load(DIR + 'latlon/GK2A_Lon_2km.npy')
GK2A_Lat=np.load(DIR + 'latlon/GK2A_Lat_2km.npy')

missing_data=np.zeros((New_IBTrACS.shape[0], 4))
ct=0
for L in range(New_IBTrACS.shape[0]):
#for L in np.arange(3290, 3806):        
    Year=str(int(New_IBTrACS[L, 2]))
    Month=str(int(New_IBTrACS[L, 3])).zfill(2)
    Day=str(int(New_IBTrACS[L, 4])).zfill(2)
    Hour=str(int(New_IBTrACS[L, 5])).zfill(2)
    
    Date=Year + Month + Day
    
    Out_DIR=DIR + 'L1B/TC_Area/' + Year + '/' + Month
    os.makedirs(Out_DIR, exist_ok=True)


    IR105_path=DIR + 'L1B/FD/IR105/' + Year + Month + '/' + Day + '/' + Hour + '/gk2a_ami_le1b_ir105_fd020ge_' + Date + Hour + '00.nc'
    IR123_path=DIR + 'L1B/FD/IR123/' + Year + Month + '/' + Day + '/' + Hour + '/gk2a_ami_le1b_ir123_fd020ge_' + Date + Hour + '00.nc'
    SW038_path=DIR + 'L1B/FD/SW038/' + Year + Month + '/' + Day + '/' + Hour + '/gk2a_ami_le1b_sw038_fd020ge_' + Date + Hour + '00.nc'
    WV069_path=DIR + 'L1B/FD/WV069/' + Year + Month + '/' + Day + '/' + Hour + '/gk2a_ami_le1b_wv069_fd020ge_' + Date + Hour + '00.nc'
    
    if (Hour=='06') & (os.path.exists(IR105_path)==False):

        IR105_path=DIR + 'L1B/FD/IR105/' + Year + Month + '/' + Day + '/' + Hour + '/gk2a_ami_le1b_ir105_fd020ge_' + Date + Hour + '10.nc'
        IR123_path=DIR + 'L1B/FD/IR123/' + Year + Month + '/' + Day + '/' + Hour + '/gk2a_ami_le1b_ir123_fd020ge_' + Date + Hour + '10.nc'
        SW038_path=DIR + 'L1B/FD/SW038/' + Year + Month + '/' + Day + '/' + Hour + '/gk2a_ami_le1b_sw038_fd020ge_' + Date + Hour + '10.nc'
        WV069_path=DIR + 'L1B/FD/WV069/' + Year + Month + '/' + Day + '/' + Hour + '/gk2a_ami_le1b_wv069_fd020ge_' + Date + Hour + '10.nc'
             
    
    if os.path.exists(IR105_path)==True:
        
        IR105_file=Dataset(IR105_path, 'r', format='netcdf4')
        IR123_file=Dataset(IR123_path, 'r', format='netcdf4')
        SW038_file=Dataset(SW038_path, 'r', format='netcdf4')
        WV069_file=Dataset(WV069_path, 'r', format='netcdf4')

        # %% Cut select TC Area
        TC_Lon=New_IBTrACS[L, 6]
        TC_Lat=New_IBTrACS[L, 7]
        
        diff=np.sqrt( ((GK2A_Lon-TC_Lon)**2) + ((GK2A_Lat-TC_Lat)**2) )
        
        x=np.where(diff==np.nanmin(diff))[0][0]
        y=np.where(diff==np.nanmin(diff))[1][0]
        
        cut_lat=GK2A_Lat[x-Pixels:x+Pixels+1, y-Pixels:y+Pixels+1]
        cut_lon=GK2A_Lon[x-Pixels:x+Pixels+1, y-Pixels:y+Pixels+1]
    
        IR105=np.array(IR105_file.variables['image_pixel_values'][x-Pixels:x+Pixels+1, y-Pixels:y+Pixels+1], dtype=float)
        IR123=np.array(IR123_file.variables['image_pixel_values'][x-Pixels:x+Pixels+1, y-Pixels:y+Pixels+1], dtype=float)
        SW038=np.array(SW038_file.variables['image_pixel_values'][x-Pixels:x+Pixels+1, y-Pixels:y+Pixels+1], dtype=float)
        WV069=np.array(WV069_file.variables['image_pixel_values'][x-Pixels:x+Pixels+1, y-Pixels:y+Pixels+1], dtype=float)
        
        IR105[np.where(IR105==32768)]=np.NaN
        IR123[np.where(IR123==32768)]=np.NaN
        SW038[np.where(SW038==32768)]=np.NaN
        WV069[np.where(WV069==32768)]=np.NaN
                
        del(IR105_file, IR123_file, SW038_file, WV069_file)
        del(diff, x, y, TC_Lon, TC_Lat)

# %%
        
        GK2A_table=pd.read_csv(DIR + 'L1B/manual/20191115_gk2a_ami_calibration_table_v3.1_ir133_srf_shift.csv')
        
        Digital_Count=np.array(GK2A_table['Digital_Count'])
        
        IR105_table=np.array(GK2A_table['IR105_Brightness_Temperature'])
        IR123_table=np.array(GK2A_table['IR123_Brightness_Temperature'])
        SW038_table=np.array(GK2A_table['IR038_Brightness_Temperature'])
        WV069_table=np.array(GK2A_table['IR069_Brightness_Temperature'])
        
        for D in range(IR105_table.shape[0]):
            IR105[np.where(IR105==D)]=IR105_table[D]
            IR123[np.where(IR123==D)]=IR123_table[D]
            SW038[np.where(SW038==D)]=SW038_table[D]
            WV069[np.where(WV069==D)]=WV069_table[D]
        del(D)
        
        """
        Nor_IR105=IR105/np.where(np.isnan(IR105_table)==0)[0][-1]
        Nor_IR123=IR123/np.where(np.isnan(IR123_table)==0)[0][-1]
        Nor_SW038=SW038/np.where(np.isnan(SW038_table)==0)[0][-1]
        Nor_WV069=WV069/np.where(np.isnan(WV069_table)==0)[0][-1]
        """
        del(GK2A_table, Digital_Count)
        del(IR105_table, IR123_table, SW038_table, WV069_table)


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
        TC_info[3:, 1]=New_IBTrACS[L, 6:28]

        
        NCFILE=Dataset(Out_DIR + '/GK2A_AMI_LE1B_TC_AREA_020GE_rot_000_' + Date + Hour + '00_' + str(int(New_IBTrACS[L, 0])) + '.nc',
                       'w', format='NETCDF4')
        
        info_i=NCFILE.createDimension('info_i', TC_info.shape[0])
        info_j=NCFILE.createDimension('info_j', TC_info.shape[1])
        lat=NCFILE.createDimension('Latitude', cut_lat.shape[0])
        lon=NCFILE.createDimension('Longitude', cut_lon.shape[0])
        
        lats=NCFILE.createVariable('lat', 'f8', ('Longitude', 'Latitude'))
        lons=NCFILE.createVariable('lon', 'f8', ('Longitude', 'Latitude'))

        ir105=NCFILE.createVariable('IR105', 'f8', ('Longitude', 'Latitude'))
        ir123=NCFILE.createVariable('IR123', 'f8', ('Longitude', 'Latitude'))
        sw038=NCFILE.createVariable('SW038', 'f8', ('Longitude', 'Latitude'))
        wv069=NCFILE.createVariable('WV069', 'f8', ('Longitude', 'Latitude'))

        #nor_ir105=NCFILE.createVariable('Nor_IR105', 'f8', ('Longitude', 'Latitude'))
        #nor_ir123=NCFILE.createVariable('Nor_IR123', 'f8', ('Longitude', 'Latitude'))
        #nor_sw038=NCFILE.createVariable('Nor_SW038', 'f8', ('Longitude', 'Latitude'))
        #nor_wv069=NCFILE.createVariable('Nor_WV069', 'f8', ('Longitude', 'Latitude'))
        
        tc_info=NCFILE.createVariable('TC_info', 'S8', ('info_i', 'info_j'))

        NCFILE.description = 'Extract only TC area from GK2A AMI LE1B 2km data'
        NCFILE.history = 'Created on ' + TT.ctime(TT.time())
        NCFILE.author = 'Y-H Baek, Typhoon research center, Jeju national university'
        lats.units = 'degrees'
        lons.units = 'degrees'
        ir105.units = 'IR105 Brightness temperature (TBB), unit:K'
        ir123.units = 'IR123 Brightness temperature (TBB), unit:K'
        sw038.units = 'SW038 Brightness temperature (TBB), unit:K'
        wv069.units = 'WV069 Brightness temperature (TBB), unit:K'

        #nor_ir105.units = 'Normalize by dividing by the largest of the IR105 correction values(8152)'
        #nor_ir123.units = 'Normalize by dividing by the largest of the IR123 correction values(8154)'
        #nor_sw038.units = 'Normalize by dividing by the largest of the SW038 correction values(16344)'
        #nor_wv069.units = 'Normalize by dividing by the largest of the WV069 correction values(8152)'
        
        tc_info.units = 'MWS : kt, DIST2LAND : km, Pres : hPa, POCI : hPa, EYE : n mi, RMW : n mi, ROCI : n mi, R34 : n mi, R50 : n mi, R64 : n mi, Heading angle : degrees'

        lats[:]=cut_lat
        lons[:]=cut_lon
        
        ir105[:]=IR105
        ir123[:]=IR123
        sw038[:]=SW038
        wv069[:]=WV069
        
        #nor_ir105[:]=Nor_IR105
        #nor_ir123[:]=Nor_IR123
        #nor_sw038[:]=Nor_SW038
        #nor_wv069[:]=Nor_WV069
        
        tc_info[:]=TC_info
        
        NCFILE.close()
        del(TC_info, info_i, info_j, lat, lon)
        del(lats, lons, ir105, ir123, sw038, wv069)
        #del(nor_ir105, nor_ir123, nor_sw038, nor_wv069)
        del(tc_info)
        del(cut_lat, cut_lon)
        del(IR105, IR123, SW038, WV069)
        #del(Nor_IR105, Nor_IR123, Nor_SW038, Nor_WV069)
    # %%
    else:
        missing_data[ct, 0]=Year
        missing_data[ct, 1]=Month
        missing_data[ct, 2]=Day
        missing_data[ct, 3]=Hour
        ct=ct+1

    del(IR105_path, IR123_path, SW038_path, WV069_path)

    print(L, '/', str(New_IBTrACS.shape[0]), Date, '-', Hour)
    
    del(Year, Month, Day, Hour, Date)

del(L)
del(GK2A_Lon, GK2A_Lat, Pixels, New_IBTrACS, ct)
 
missing_data=missing_data[np.where(missing_data[:, 0]!=0)[0], :]
np.save(DIR + 'L1B/missing_data_GK2A_TCvitals.npy', missing_data)
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




