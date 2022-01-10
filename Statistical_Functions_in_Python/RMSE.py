#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 13:43:10 2021

@author: yhbaek
"""
import numpy as np

def rmse(predictions, targets): 

    return np.sqrt(((predictions - targets) ** 2).mean())
