# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 23:06:02 2019

"""
import math

def leftMultVecMat(vec, Mat):
    """
    Inputs: vector, Matrix
    Outputs: left side multiplication (vec*Mat)
    """
    result = []
    for i in range(len(Mat[0])): #this loops through columns of the matrix
        total = 0
        for j in range(len(vec)): #this loops through vector coordinates & rows of matrix
            total += vec[j] * Mat[j][i]
        result.append(total)
    return result

def rightMultVecMat(vec, Mat):
    """
    Inputs: vector, Matrix
    Outputs: right side multiplication (Mat*vec)
    """
    result = []
    for i in range(len(Mat[0])): #this loops through columns of the matrix
        total = 0
        for j in range(len(vec)): #this loops through vector coordinates & rows of matrix
            total += vec[j] * Mat[i][j]
        result.append(total)
    return result

def normalize(vec):
    """
    Inputs: vector
    Outputs: normalized vector
    """
    norm =  normVec(vec)
    vec = [vec[0]/norm, vec[1]/norm, vec[2]/norm]
    return vec

def normVec(vec):
    """
    Inputs: vector
    Outputs: norm of vector
    """
    return math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

def dot(v1, v2):
    """
    Inputs: 2 vectors
    Outputs: dot product of 2 vectors
    """
    return sum(x*y for x,y in zip(v1,v2))

def scaleVec(inVec,scalar):
    """
    Inputs: vector and scalar multiplier value
    Outputs: vector with each element multiplied by the input scalar
    """
    outVec = [inVec[0]*scalar, inVec[1]*scalar, inVec[2]*scalar]
    return outVec

def subtractVec(vec1,vec2):
    """
    Inputs: 2 vectors
    Outputs: vector1 - vector2 in vector form
    """
    return [vec1[0] - vec2[0], vec1[1] - vec2[1], vec1[2] - vec2[2]]

def addVec(vec1,vec2):
    """
    Inputs: 2 vectors
    Outputs: vector1 + vector2 in vector form
    """
    return [vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2]]

def sense2vector(measurements, r_Earth2Sun, r_sat):
    """
    Inputs: 
        measurements: raw measurement values from 6 sun sensors
        r_Earth2Sun: Earth to Sun vector
        r_sat: position of satellite in ECI
    Outputs:
        sat2sun: satellite to sun 3-vector
    """
    
    #unpack measurements  
    sun_sense1 = measurements[0]
    sun_sense2 = measurements[1]
    sun_sense3 = measurements[2]
    sun_sense4 = measurements[3]
    sun_sense5 = measurements[4]
    sun_sense6 = measurements[5]
    
    irrad_vec = [sun_sense1 - sun_sense2,  sun_sense3 - sun_sense4, sun_sense5 - sun_sense6] #create irradiance vector from sensor values
    irrad_vec = normalize(irrad_vec) #normalize irradiance vector
    
    #need to convert to body frame without numpy
    albedo = scaleVec(normalize(r_sat), 0.2)
    sat2sun = normalize(subtractVec(irrad_vec,albedo)) #vector subt. irradiance vec and albedo vec, normalize
    return sat2sun

def vector2sense(sat2sun, r_sat):
    irrad_vec = addVec(sat2sun,scaleVec(r_sat))
    delta1 = irrad_vec[0]
    delta2 = irrad_vec[1]
    delta3 = irrad_vec[2]
    #need to figure out minimum value and do body frame conversion
    

def isEclispe(r_sat, r_Earth2Sun, Re):
    """
    Inputs: 
        r_sat: position of satellite in ECI
        r_Earth2Sun: Earth to Sun vector
        Re: Radius of Earth
    Outputs:
        True if satellite is in eclipse with Earth
        False if satellite is NOT in eclipse with Earth
    """
    r_Sat2Sun = subtractVec(r_Earth2Sun,r_sat)
    theta = math.acos(dot(r_Sat2Sun,scaleVec(r_sat,-1))/(normVec(r_Sat2Sun)*normVec(r_sat)))
    transit_L = normVec(r_sat)*math.sin(theta)
    if transit_L > Re:
        return False
    else:
        return True
vec = [1,2,3]
Mat = [[1,2,3],[1,2,3],[1,2,3]]
print(rightMultVecMat(vec, Mat))
    
a = dot([1,2,3], [4,5,6])
b = norm([1,2,3])
print(a/b)


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





