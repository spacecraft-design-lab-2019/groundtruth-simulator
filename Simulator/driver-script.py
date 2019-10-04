# -*- coding: utf-8 -*-
"""
Driver Script

Initializes and runs groundtruth simulator.

Pseudocode:
    -initialize with starting position & velocity

    LOOP:
        -calculate environment
        -calculate forces and torques (dynamics)
        -integrate (kinematics)
        -obtain new position/velocity

        -model sensors with noise

        -FEED to controller
        -take input from controller
        -add noise
        -account for control input in dynamics
"""

#-------- Initialize Workspace and Variables ------------
import numpy as np
from initialize import *
from rk4_step import rk4_step
from calc_statedot import calc_statedot


#------------ Run Simulation -----------------------------
t = np.arange(tspan[0], tspan[1]+tstep, tstep)
state = np.zeros((np.shape(t)[0], np.shape(state_initial)[0]))

for idx in range(t.shape[0]):
    # Integrate
    ti, state[idx+1, :] = rk4_step(calc_statedot, t[idx], state[idx, :], tstep)
    
    # Normalize the Quaternion Vector
    state[idx+1, 3:7] = state[idx+1, 3:7]/np.linalg.norm(state[idx+1, 3:7])

#-------------------- Plot------------------------------------
# can add plotting/analysis scripts here