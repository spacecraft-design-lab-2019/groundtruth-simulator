# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 14:02:56 2019

@author: Kevin
"""
import numpy as np
import pandas as pd
import os

#directory
dir_path = os.path.dirname(os.path.realpath(__file__))

Re = 6378.1378366 #km, measured at equator
mSun = 1.98855e30 #kg
wEarth = 2*np.pi/86164.1 #orbital speed earth
GM = 3.986e5


#Spherical Gravity
## NOTE: the below xlsx file contains Gravity Coefficients for an asteroid, NOT Earth
gravity_coefficients = pd.read_excel("lookups/GravityCoefficients.xlsx",header=None)
order = 35
degree = 35
C = np.zeros((order,degree)) #not sure if the order and degree are in the right spots
S = np.zeros((order,degree))

for i in range(0,len(gravity_coefficients)):
    ord = gravity_coefficients.iloc[i,0]
    deg = gravity_coefficients.iloc[i,1]
    C[ord,deg] = gravity_coefficients.iloc[i,2]
    S[ord,deg] = gravity_coefficients.iloc[i,3]
