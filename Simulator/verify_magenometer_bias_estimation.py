# -*- coding: utf-8 -*-

"""
Verify detumble algorithms
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import sim_config as config
from simulator import Simulator
from mpl_toolkits.mplot3d import Axes3D

# import GNC functions
import sys

sys.path.append('/home/eleboeuf/Documents/GNC')
import detumble_cpp as dcpp


# -----------------------Clear Figures----------------------------
# clear figures
plt.close('all')

# -----------------Configuration / Parameters--------------------
tspan = np.array([0, 1])  # [sec]
L_cmd = np.zeros(3)  # initially command 0 torque

# ----------------Initialize / Setup Workspace------------------
# setup sim
sim = Simulator(config)
max_dipoles = np.array([[8.8e-3], [1.373e-2], [8.2e-3]])  # TODO: Add dipole saturation to spacecraft properties

# preallocate memory
T = np.arange(0, tspan[1] + config.tstep, config.tstep)
B_body_history = np.zeros((np.shape(T)[0], 3))

# ---------------------Propagate---------------------------
t = time.time()
for i in range(len(T)):
    # Simulator
    sensors = sim.step(config.tstep, L_cmd)

    # command torque based on sensors (currently no noise addition 11/17)
    B_sensed = sensors[0:3]
    B_body_history[i + 1, :] = np.transpose(B_sensed)

elapsed = time.time() - t
print(elapsed)
# ------------------------Plot-----------------------------


with plt.rc_context(rc={'interactive': False}):
    plt.show()


