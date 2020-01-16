# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 23:06:02 2019

"""
import math
import numpy as np
import conversions as conv
import astropy_sun_position as sun

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
    if all([v == 0 for v in vec]): #returns 0 vector if vec = [0,0,0]
        return vec
    else:
        return [x/mag for x in vec]

def transpose(M):
    I = range(len(M))
    J = range(len(M[0]))
    return [[M[i][j] for i in I] for j in J]

def norm(vec):
    return math.sqrt(dot(vec, vec))

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

def sense2vector(meas, r_sat, q_eci2body, albedo = True):
    """
    Inputs:
        meas: raw measurement values from 6 sun sensors. Arranged: [x, -x, y, -y, z, -z]
        r_Earth2Sun: Earth to Sun vector
        r_sat: position of satellite in ECI
    Outputs:
        sat2sun: satellite to sun 3-vector (in body frame)
    """

    irrad_vec = [meas[0] - meas[1],  meas[2] - meas[3], meas[4] - meas[5]] #create irradiance vector from sensor values
    irrad_vec = normalize(irrad_vec) # normalize irradiance vector

    if albedo:
        alb = conv.quatrot(q_eci2body, scale(normalize(r_sat), 0.2)) #convert to body frame
        sat2sun = normalize(sub(irrad_vec, alb)) #vector subt. irradiance vec and albedo vec, normalize
    else:
        sat2sun = irrad_vec

    return sat2sun #in body frame

def vector2sense(sat2sun, r_sat, q_eci2body, albedo = True):
    """
    Convert the true sun vector to an array of voltages corresponding to sun sensor measurements

    Inputs:
        sat2sun: satellite to sun vector in body
        r_sat: position of satellite in eci
    """
    sat2sun = normalize(sat2sun) #make sure vec is normalized
    if albedo:
        alb = conv.quatrot(q_eci2body, scale(normalize(r_sat), 0.2)) #body frame
        irrad_vec = normalize(add(sat2sun, alb) )
    else: 
        irrad_vec = sat2sun
        
    return deltas2measure(irrad_vec)

def deltas2measure(deltas):
    """
    Helper function to arrange vector of measurements from adjusted light vector
    Can create a measurement list regardless of input list length
    """
    measurements = []

    for d in deltas:
        if d < 0:
            measurements += [0, abs(d)]
        else:
            measurements += [abs(d), 0]

    return np.array(measurements)


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



#in_meas = [1, 0, 0, 0, 0, 0]
#r_sat = [8000, 0, 0]
#q = [1, 0, 0, 0] # identity quaternion
#sat2sun = sense2vector(in_meas, r_sat, q, albedo = True)
#out_meas = vector2sense(sat2sun, r_sat, conv.conj(q), albedo = True)
#print(out_meas)
#
#meas = [1, 2, 0, 0, 0, 0]
#r_sat = [8000, 0, 0]
#q = [1, 0, 0, 0] # identity quaternion
#
#vec = sense2vector(meas, r_sat, q, albedo = False)
#print(vec)




