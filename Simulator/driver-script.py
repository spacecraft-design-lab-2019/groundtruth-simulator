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
from simulator import *
import matplotlib.pyplot as plt


#------------ Run Simulation -----------------------------
t = np.arange(tspan[0], tspan[1]+tstep, tstep)
state = np.zeros((np.shape(t)[0], np.shape(state_initial)[0]))
state[0, :] = state_initial

for idx in range(t.shape[0]-1):
    # Integrate
    ti, state[idx+1, :] = rk4_step(calc_statedot, t[idx], state[idx, :], tstep)
    
    # Normalize the Quaternion Vector
    state[idx+1, 3:7] = state[idx+1, 3:7]/np.linalg.norm(state[idx+1, 3:7])

#-------------------- Plot------------------------------------
# can add plotting/analysis scripts here
plt.plot(state[:, 0], state[:, 1])
plt.show()