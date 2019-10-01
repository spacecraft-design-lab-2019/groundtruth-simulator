# -*- coding: utf-8 -*-
"""
Initialization Script

Sets simulation parameters and initializes spacecraft. Run by driver-script
before beginning simulation.

NOTE: based on \Reference\Hridu Old Simulator\Model Scripts\initialize.m
"""
#Importing libraries
import numpy as np

# Simulation Parameters
class SpacecraftParams():
    """
    A class to store all spacecraft parameters
    """
    def __init__(self,
                I = np.array([[17,0,0],[0,18,0],[0,0,22]]),
                r = np.array([1,1,1]),
                v = np.array([1,1,1])):
        self.I = I
        self.r = r
        self.v = v

#testing SpacecraftParams
s = SpacecraftParams()
print(s.I)
s.r = np.array([5,5,5])
print(s.r)
s.v[0] = 9
print(s.v)
# Initial Spacecraft State Vector


# Mass Properties of Spacecraft


# Structure of Spacecraft


# Magnetorquers


# Sensors


# EKF Parameters


# Global Constants
