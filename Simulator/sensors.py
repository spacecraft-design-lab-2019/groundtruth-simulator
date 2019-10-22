# -*- coding: utf-8 -*-
"""
Sensor modeling for satellite flight simulator
incoprorates a scale factor, axial misalignment, slow varying and fast varying offsets
"""

import numpy as np
pi = (np.pi)/6

def linear_error_model(x, mean, cov, mis_angle, B, V):
    y = scale_factor(mean,cov)*misalign(mis_angles)*x + B + V
    return y

def scale_factor(mean,cov):
    scaleF = np.eye(3) + np.diag(np.random.multivariate_normal(mean,cov,1)[0])
    return scaleF

def misalign(mis_angles):
    M12 = mis_angles[0]
    M13 = mis_angles[1]
    M23 = mis_angles[2]
    M21 = M12
    M32 = M23
    M31 = M13
    Mis_Mat = np.array([[1, M12, M13],[M21, 1, M23],[M31, M32, 1]])
    return Mis_Mat
    
    
mean = np.array([0, 0, 0])
cov = np.diag([.01, .01, .01])
mis_angles = np.random.multivariate_normal([pi/6,pi/6,pi/6],np.diag([pi/10,pi/10,pi/10]),1)[0]
scaleF = linear_error_model(mean,cov)
print(scaleF)