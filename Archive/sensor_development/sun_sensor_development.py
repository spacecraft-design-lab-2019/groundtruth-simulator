# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:26:23 2020

@author: Kevin
archived code for sensor development
"""

"""
Must have sun_sensor_math file
Test Code Under here:
"""
from sun_sensor_math import *
deltas = [1,-1,1]
print(deltas2measure(deltas))
print(deltas2measure(deltas) == [1,0,0,1,1,0])

vec = [1.09, -2.02, 3.04]
scalar = 0
print(scale(vec, scalar))

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
