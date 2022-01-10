#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 22:23:07 2021

@author: yhbaek
"""

import os
import sys
import datetime
import numpy as np

sys.path.append('/data/GK2A/')
from Function_of_API_gk2a_download import down_gk2a

# COMS : VIS : 0.67, SWIR : 2.7, WV : 6.7, IR1 : 10.8, IR2 : 12.0
lv='LE1B'
AREA='FD'
#CH='WV069'	# SW038, WV069, IR105, IR123
DIR='/data/GK2A/'



# IOT : 2019.07.01-2020.03.31
# %%
# missing : 2020, 2, 26 6~23H
# missing : 2020, 3, 28 6~8H
# 2021, 10, 25
Str_date=datetime.date.toordinal(datetime.date(2021, 10, 26))
End_date=datetime.date.toordinal(datetime.date(2021, 12, 31))


for L in np.arange(Str_date, End_date+1, 1):
    date=datetime.date.fromordinal(L)
    for Hour in np.arange(0, 24, 1):
        if Hour!=6:
            time=datetime.date.strftime(date, "%Y%m%d") + str(Hour).zfill(2) + '00'
        elif Hour==6:
            time=datetime.date.strftime(date, "%Y%m%d") + str(Hour).zfill(2) + '10'
			
        down_gk2a(time, lv, 'SW038', AREA, DIR)
        down_gk2a(time, lv, 'WV069', AREA, DIR)   
        down_gk2a(time, lv, 'IR105', AREA, DIR)
        down_gk2a(time, lv, 'IR123', AREA, DIR)

        del(time)
    del(Hour)
    del(date)
del(L)

del(Str_date, End_date)    

"""

D=datetime.date.toordinal(datetime.date(2019, 12, 4))
date=datetime.date.fromordinal(D)
time=datetime.date.strftime(date, "%Y%m%d") + '1400'

#down_gk2a(time, lv, 'SW038', AREA, DIR)
#down_gk2a(time, lv, 'WV069', AREA, DIR)
down_gk2a(time, lv, 'IR105', AREA, DIR)
down_gk2a(time, lv, 'IR123', AREA, DIR)

"""







