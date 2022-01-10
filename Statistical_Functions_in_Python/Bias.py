#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 00:45:17 2021

@author: yhbaek
"""

import numpy as np

def bias(predictions, targets): 

    return (predictions - targets).mean()