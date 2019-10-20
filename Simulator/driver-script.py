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
import matplotlib.pyplot as plt


#------------------ Run Simulation -----------------------------

magnes = SpacecraftState()
world = Environment()

t = np.arange(magnes.t, tspan[1]+tstep, tstep)
state_history = np.zeros((np.shape(t)[0], np.shape(magnes.state())[0]))
state_history[0, :] = magnes.state()

for idx in range(t.shape[0]-1):
    # Integrate
    magnes = rk4_step(calc_statedot, magnes, world, tstep)

    # Normalize the Quaternion Vector
    magnes.normalize_quat();

    # Save to History
    state_history[idx+1, :] = magnes.state()


#-------------------- Plot------------------------------------

plt.plot(state_history[:, 0], state_history[:, 1])
plt.xlabel("X_ECI")
plt.ylabel("Y_ECI")
plt.grid()
plt.show()
