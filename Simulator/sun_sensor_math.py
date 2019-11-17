# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 23:06:02 2019

"""
import math

def normalize(vec):
    norm =  normVec(vec)
    vec = [vec[0]/norm, vec[1]/norm, vec[2]/norm]
    return vec

def normVec(vec):
    return math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

def dot(v1, v2):
    return sum(x*y for x,y in zip(v1,v2))

def scaleVec(inVec,scaler):
    outVec = [inVec[0]*scaler, inVec[1]*scaler, inVec[2]*scaler]
    return outVec

def subtractVec(vec1,vec2):
    return [vec1[0] - vec2[0], vec1[1] - vec2[1], vec1[2] - vec2[2]]

def sense2vector(measurements, r_Earth2Sun, r_sat):
    sun_sense1 = measurements[0]
    sun_sense2 = measurements[1]
    sun_sense3 = measurements[2]
    sun_sense4 = measurements[3]
    sun_sense5 = measurements[4]
    sun_sense6 = measurements[5]
    
    irrad_vec = [sun_sense1 - sun_sense2,  sun_sense3 - sun_sense4, sun_sense5 - sun_sense6]
    irrad_vec = normalize(irrad_vec)
    
    albedo = scaleVec(normalize(r_Earth2Sun), 0.2)
    sat2sun = irrad_vec
    return sat2sun

def isEclispe(r_sat, r_Earth2Sun, Re):
    r_Sat2Sun = subtractVec(r_Earth2Sun,r_sat)
    theta = math.acos(dot(r_Sat2Sun,scaleVec(r_sat,-1))/(normVec(r_Sat2Sun)*normVec(r_sat)))
    transit_L = normVec(r_sat)*math.sin(theta)
    if transit_L > Re:
        return True
    else:
        return False
    
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





