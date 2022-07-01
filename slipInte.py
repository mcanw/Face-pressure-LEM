# -*- coding: utf-8 -*-
"""
Created on Thu May 27 09:14:59 2021

@author: hdec
"""
import numpy as np
import math
def interFric(fric):
    #calculate the cofficient of moblized reduction coefficient at tunnel interface
    radiansfric = math.radians(fric)
    radiansdelta = math.atan(math.sin(radiansfric)*math.cos(radiansfric)/(1+math.sin(radiansfric)**2))
    inteRed = math.tan(radiansdelta)/math.tan(radiansfric)
    # print(inteRed,radiansfric,radiansdelta)
    return(inteRed)