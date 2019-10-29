# -*- coding: utf-8 -*-
"""
Initialization Script

Sets simulation parameters and initializes spacecraft. Run by driver-script
before beginning simulation.

"""

#-------------------------Setup---------------------------------

# Importing libraries
import numpy as np

# Simulation Parameters
tspan = np.array([0, 864])    # [sec]
tstep = .1                     # [sec] - 10 Hz
mjd_start = 58777.740671


# Initial Spacecraft State
r_i = np.array([6712880.93e-3,1038555.54e-3,-132667.04e-3])
q_i = np.array([1, 0, 0, 0])
v_i = np.array([-831.937369e-3,4688.525767e-3,-6004.570270e-3])
w_i = np.array([.1, .5, -.3])
state_i = np.r_[r_i, q_i, v_i, w_i]


#--------------------------------------------------------------