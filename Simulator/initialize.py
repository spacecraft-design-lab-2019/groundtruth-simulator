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
tspan = np.array([0, 5])    # [sec]
tstep = .1                  # [sec] - 10 Hz


# Initial Spacecraft State Vector
# ------ can later use oe2eci to get starting r and v based on injection orbit
r_initial = np.array([0,1.5e11,0])
v_initial = np.array([30000,0,0])
q_initial = np.array([0, 0, 0, 1])
w_initial = np.array([0, 0, 0])
state_initial = np.r_[r_initial, q_initial, v_initial, w_initial]


# Mass Properties of Spacecraft
class SpacecraftParams():
    """
    A class to store all spacecraft parameters
    """
    def __init__(self,
                I = np.array([[17,0,0],[0,18,0],[0,0,22]])):
        self.I = I
        

# Structure of Spacecraft


# Magnetorquers


# Sensors


# EKF Parameters


# Global Constants
mu = 6.67408e-11*2e30

