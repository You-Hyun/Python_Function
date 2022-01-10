#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 22:23:07 2021

@author: yhbaek
"""

import sys
import os
import datetime
import numpy as np

sys.path.append('/home/yhbaek/STUDY/Python_function/')
from Function_of_API_gk2a_download import down_gk2a

lv='LE1B'
AREA='FD'
CH='SW038'
DIR='/home/yhbaek/'

# %%
Str_date=datetime.date.toordinal(datetime.date(2020, 5, 11))
End_date=datetime.date.toordinal(datetime.date(2020, 5, 16))

for L in np.arange(Str_date, End_date+1, 1):
    date=datetime.date.fromordinal(L)
    for Hour in np.arange(0, 24, 1):
        time=datetime.date.strftime(date, "%Y%m%d") + str(Hour).zfill(2) + '00'
    
        down_gk2a(time, lv, CH, AREA, DIR)
        del(time)
    del(Hour)
    del(date)
del(L)

del(Str_date, End_date)    

