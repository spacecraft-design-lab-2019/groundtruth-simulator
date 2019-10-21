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

#------------ Initialize Workspace and Variables ------------
import numpy as np
from initialize import *
from propagate_step import *
from constants import *
import matplotlib.pyplot as plt


#------------------ Run Simulation -----------------------------

world = Environment()
struct = SpacecraftStructure()

t = np.arange(0, tspan[1]+tstep, tstep)
state_history = np.zeros((np.shape(t)[0], np.shape(state_i)[0]))
state_history[0, :] = state_i

for idx in range(t.shape[0]-1):
    # Integrate
    state_history[idx+1, :] = rk4_step(lambda t, state: calc_statedot(t, state, world, struct), t[idx], state_history[idx, :], tstep)

    # Normalize the Quaternion Vector
    state_history[idx+1, 3:7] = state_history[idx+1, 3:7]/np.linalg.norm(state_history[idx+1, 3:7])


#-------------------- Plot------------------------------------

plt.figure()
plt.plot(state_history[:, 0], state_history[:, 1])
plt.xlabel("X_ECI")
plt.ylabel("Y_ECI")
plt.grid()
plt.show()

plt.figure()
plt.plot(t, state_history[:, 3:7])
plt.xlabel('time')
plt.ylabel('quaternions')
plt.grid()
plt.show()