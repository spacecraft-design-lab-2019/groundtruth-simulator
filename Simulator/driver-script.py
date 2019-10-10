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
from propagate_step import *
import matplotlib.pyplot as plt


#------------ Run Simulation -----------------------------

ISS = SpacecraftState()

t = np.arange(ISS.t, tspan[1]+tstep, tstep)
state_history = np.zeros((np.shape(t)[0], np.shape(ISS.state())[0]))
state_history[0, :] = ISS.state()

for idx in range(t.shape[0]-1):
    # Integrate
    ISS = rk4_step(calc_statedot, ISS, tstep)
    state_history[idx+1, :] = ISS.state()

    # Normalize the Quaternion Vector
    state_history[idx+1, 3:7] = state_history[idx+1, 3:7]/np.linalg.norm(state_history[idx+1, 3:7])

#-------------------- Plot------------------------------------
# can add plotting/analysis scripts here
plt.plot(state_history[:, 0], state_history[:, 1])
plt.xlabel("X_ECI")
plt.ylabel("Y_ECI")
plt.grid()
plt.show()
