#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 00:40:34 2021

@author: yhbaek
"""

import numpy as np

def mae(predictions, targets): 

    return (abs(predictions - targets)).mean()