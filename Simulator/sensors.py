# -*- coding: utf-8 -*-
"""
Integration of Sensor Model
"""

import numpy as np
from numpy import pi
from sensormodels import *


def getTmatrix(scaleF, caSense):
    scaleFmat = np.eye(3) + np.diag(whitenoise(scaleF,3))
    misalign = np.reshape(np.random.multivariate_normal(np.zeros(9),caSense*np.eye(9)),(3,3))
    np.fill_diagonal(misalign, 0)
    T = np.dot(scaleFmat,misalign)
    return T
    
#Define Gyro
scaleF = 0.0002
caSense = 0.02
T =  getTmatrix(scaleF, caSense) 
gyro = Sensor(errormodel = LinearErrorModel.withDim(3, b = np.ones(3), cov = 0.0005))
X = gyro.measure(np.ones(3))
print(X)
    

