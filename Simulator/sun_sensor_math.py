# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 23:06:02 2019

"""
import math
import numpy as np

#placeholder function, need to fix
def eci2body(vec, R):
    """
    Inputs: vector in eci, Rotation matrix for eci to body
    Outputs: vector in body frame
    """
    return matTimesVec(R, vec)

def vecTimesMat(x, M):
    """
    Vector multiplication for lists.
    Performs x^T * M
    """
    return [dot(m, x) for m in transpose(M)]

def matTimesVec(M, x):
    """
    Vector multiplication for lists.
    Performs M * x
    """
    return [dot(m, x) for m in M]

def normalize(vec):
    """
    Inputs: vector
    Outputs: normalized vector
    """
    mag = norm(vec)
    return [x/mag for x in vec]

def transpose(M):
    I = range(len(M))
    J = range(len(M[0]))
    return [[M[i][j] for i in I] for j in J]

def norm(vec):
    return dot(vec, vec)

def dot(v1, v2):
    return sum(x*y for x,y in zip(v1,v2))

def scale(vec, scalar):
    return [x*scalar for x in vec]

def sub(vec1, vec2):
    """
    Inputs: 2 vectors
    Outputs: vector1 - vector2 in vector form
    """
    return [x - y for x, y in zip(vec1, vec2)]

def add(vec1, vec2):
    """
    Inputs: 2 vectors
    Outputs: vector1 + vector2 in vector form
    """
    return [x + y for x, y in zip(vec1, vec2)]

def sense2vector(measurements, r_Earth2Sun, r_sat):
    """
    Inputs:
        measurements: raw measurement values from 6 sun sensors
        r_Earth2Sun: Earth to Sun vector
        r_sat: position of satellite in ECI
    Outputs:
        sat2sun: satellite to sun 3-vector
    NOTE: This function can use numpy
    """
    
    #unpack measurements  
    sun_sense1 = measurements[0] # + x
    sun_sense2 = measurements[1] # - x
    sun_sense3 = measurements[2] # + y
    sun_sense4 = measurements[3] # - y
    sun_sense5 = measurements[4] # + z
    sun_sense6 = measurements[5] # - z

    irrad_vec = [sun_sense1 - sun_sense2,  sun_sense3 - sun_sense4, sun_sense5 - sun_sense6] #create irradiance vector from sensor values
    irrad_vec = normalize(irrad_vec) #normalize irradiance vector

    albedo = eci2body(scale(normalize(r_sat), 0.2)) #convert to body frame
    sat2sun = normalize(sub(irrad_vec, albedo)) #vector subt. irradiance vec and albedo vec, normalize
    return sat2sun #in body frame

def vector2sense(sat2sun, r_sat, R_eci2body):
    """
    Convert the true sun vector to an array of voltages corresponding to sun sensor measurements

    Inputs:
        sat2sun: satellite to sun vector
        r_sat: position of satellite in eci
        R_eci2body: rotation matrix to convert to body frame
    """

    r_sat_body = dot(R_eci2body, scale(r_sat, 0.2))
    irrad_vec = add(sat2sun, r_sat_body)
    delta1 = irrad_vec[0]
    delta2 = irrad_vec[1]
    delta3 = irrad_vec[2]

    deltas = np.array([delta1,delta2,delta3])
    measurements = np.array([])

    for i in range(len(deltas)):
        if deltas[i] <0:
            measurements = np.append(measurements, 0)
            measurements = np.append(measurements, abs(deltas[i]))
        else:
            measurements = np.append(measurements, abs(deltas[i]))
            measurements = np.append(measurements, 0)

    return measurements

def isEclipse(r_sat, r_Earth2Sun, Re):
    """
    Outputs True if satellite is in eclipse with Earth, False otherwise.
    Inputs:
        r_sat: position of satellite in ECI
        r_Earth2Sun: Earth to Sun vector
        Re: Radius of Earth
    """
    r_Sat2Sun = sub(r_Earth2Sun, r_sat)
    theta = math.acos(dot(r_Sat2Sun, scale(r_sat,-1)) / (norm(r_Sat2Sun)*norm(r_sat)))
    transit_L = norm(r_sat)*math.sin(theta)
    if transit_L > Re:
        return False
    else:
        return True

def isEclipse2(measurements, thresh):
    """
    Test measurements determines if sat is in eclipse
    returns True if there is eclipse (all values are below a threshold)
    returns False if NOT in eclipse with earth
    Inputs: 
        measurements: list of values from sun sensors
        thresh: all sun sensor values need to be below this to be in eclipse
    Outputs:
        True if satellite is in eclipse with Earth
        False if satellite is NOT in eclipse with Earth
    """
    for x in measurements:
        if thresh>=x:
            return False
    return True
"""
Test Code Under here:
"""
a = isEclipse2([6,5,6,5,6,7],4)
print(a)
"""

deltas = np.array([200,-190,300])
measurements = np.array([])
for i in range(len(deltas)):
    if deltas[i] <0:
        measurements = np.append(measurements,0)
        measurements = np.append(measurements,abs(deltas[i]))
    else:
        measurements = np.append(measurements,abs(deltas[i]))
        measurements = np.append(measurements,0)

print(measurements)
vec = [1,2,3]
Mat = [[1,2,3],[1,2,3],[1,2,3]]
print(rightMultVecMat(vec, Mat))

a = dot([1,2,3], [4,5,6])
b = normalize([1,2,3])



sun_sense1 = 123
sun_sense2 = 231
sun_sense3 = 421
sun_sense4 = 546
sun_sense5 = 345
sun_sense6 = 112

gain = 1
sensors = [sun_sense1, sun_sense2, sun_sense3, sun_sense4, sun_sense5, sun_sense6];
for i in range(len(sensors)):
    sensors[i] = gain*sensors[i]

vector1 = [sun_sense1 - sun_sense2,  sun_sense3 - sun_sense4, sun_sense5 - sun_sense6]
vector2 = normalize(vector1)

"""



