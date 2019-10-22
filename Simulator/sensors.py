# -*- coding: utf-8 -*-
"""
Sensor modeling for satellite flight simulator
incoprorates a scale factor, axial misalignment, slow varying and fast varying offsets
"""

import numpy as np
from numpy import pi

def lin_error_mod(x, mean, cov, mis_angle, B, V):
    scaleF = scale_factor(mean,cov)
    misMat = misalign(mis_angle)
    y = np.matmul(scaleF,np.matmul(misMat,x)) + B + V
    return y

def scale_factor(mean,cov):
    scaleF = np.eye(3) + np.diag(np.random.multivariate_normal(mean,cov,1)[0])
    return scaleF

def misalign(mis_angles):
    M12 = mis_angles[0] #will need to be means and covariances for this
    M13 = mis_angles[1]
    M23 = mis_angles[2]
    M21 = M12
    M32 = M23
    M31 = M13
    Mis_Mat = np.array([[1, M12, M13],[M21, 1, M23],[M31, M32, 1]])
    return Mis_Mat
    
    
mean = np.array([0, 0, 0])
cov = np.diag([.01, .01, .01])
mis_angle = np.random.multivariate_normal([pi/6,pi/6,pi/6],np.diag([pi/10,pi/10,pi/10]),1)[0]
B = np.zeros([3,1])
V = np.zeros([3,1])
x = np.array([[1],[2],[3]])

y = lin_error_mod(x,mean,cov,mis_angle,B,V)
print(y)