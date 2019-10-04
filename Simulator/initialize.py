# -*- coding: utf-8 -*-
"""
Initialization Script

Sets simulation parameters and initializes spacecraft. Run by driver-script
before beginning simulation.

NOTE: based on \Reference\Hridu Old Simulator\Model Scripts\initialize.m
"""

# Importing libraries
import numpy as np

# Simulation Parameters
tspan = np.array([0, 8640])    # [sec]
tstep = .1                  # [sec] - 10 Hz

# Mass Properties of Spacecraft
class SpacecraftParams():
    """
    A class to store all spacecraft parameters
    """
    def __init__(self,
                I = np.array([[17,0,0],[0,18,0],[0,0,22]]),
                r = np.array([6712880.93e-3,1038555.54e-3,-132667.04e-3]),
                v = np.array([-831.937369e-3,4688.525767e-3,-6004.570270e-3]),
                q = np.array([1, 0, 0, 0]), # scalar first
                w = np.array([0, 0, 0])):
        self.I = I
        self.r = r
        self.v = v
        self.q = q
        self.w = w

# Initial Spacecraft State Vector
# ------ can later use oe2eci to get starting r and v based on injection orbit

ISS = SpacecraftParams()
state_initial = np.r_[ISS.r, ISS.q, ISS.v, ISS.w]


# Structure of Spacecraft


# Magnetorquers


# Sensors


# EKF Parameters


# Global Constants
