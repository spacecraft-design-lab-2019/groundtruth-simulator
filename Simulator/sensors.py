# -*- coding: utf-8 -*-
"""
Integration of Sensor Model
"""

import numpy as np
from numpy import pi
from sensormodels import *
from sun_model import *
from numpy import linalg as LA


def getTmatrix(scaleF, caSense):
    """
    Inputs:
        scaleF: scale factor for a particular sensor
        caSense: cross-axis sensitivity for particular sensor
    Outputs:
        T-matrix, combined misalignment and scaling matrix for linear error model
    """
    scaleFmat = np.eye(3) + np.diag(whitenoise(scaleF,3))
    misalign = np.reshape(np.random.multivariate_normal(np.zeros(9),caSense*np.eye(9)),(3,3))
    np.fill_diagonal(misalign, 0)
    T = np.dot(scaleFmat,misalign)
    return T

def sunVector(rsat,jdate):
    """
    Inputs:
        rsat: ECI position of satellite
        jdate: mean julian date time
    Output:
        position vector (3-vector) from satellite to sun
    """
    earth2sun = sun1(jdate)
    rsun = (earth2sun - rsat)/LA.norm(earth2sun - rsat)
    return rsun

#Gyro Function
def gyroModel(X, scaleF = 0.002, caSense = 0.02):
    """
    Inputs: 
        X: 3-vector for position/attitude
        scaleF: Scale Factor value for gyro
        caSense: Cross-axis sensitivity value for gyro
    Output:
        measurement model in 3-vector form
    """
    Tmat =  getTmatrix(scaleF, caSense) 
    gyro = Sensor(errormodel = LinearErrorModel.withDim(3, T = Tmat, b = np.ones(3), cov = 0.0005))
    return gyro.measure(X)

#Sun Sensor Function
def sunSenseModel(X, jdate, scaleF = 0, caSense = 0.02,):
    """
    Inputs: 
        X: 3-vector for position/attitude
        scaleF: Scale Factor value for sun sensor
        caSense: Cross-axis sensitivity value for sun sensor
    Output:
        measurement model in 3-vector form
    """
    Tmat =  getTmatrix(scaleF, caSense) 
    sunSense = Sensor(errormodel = LinearErrorModel.withDim(3, T = Tmat, b = np.ones(3), cov = 0.0005))
    return sunSense.measure(X)


#Test Code    
X = np.ones(3)
print(gyroModel(X, scaleF = 0.002, caSense = 0.02))
#print(X)

print(sunVector(np.ones(3),1929))



