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
import sys, os
dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, dir+'/GNC/')
import detumble_cpp as dcpp


# -----------------------Clear Figures----------------------------
# clear figures
plt.close('all')

# -----------------Configuration / Parameters--------------------
tspan = np.array([0, 1200])  # [sec]
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
    sensors = sim.step(L_cmd)

    # command torque based on sensors (currently no noise addition 11/17)
    B_sensed = sensors[0:3]
    B_body_history[i, :] = np.transpose(B_sensed)

elapsed = time.time() - t
print(elapsed)
#-------------------------Estimate magnetometer bias-------
bias_est = dcpp.get_bias_estimate(B_body_history)

print(sim.sensors.magnetometer.errormodel.bias_current) # true bias
print(bias_est)

# ------------------------Plot-----------------------------
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(B_body_history[:,0], B_body_history[:,1], B_body_history[:,2], 'bo')
ax.plot3D([0, bias_est[0]], [0, bias_est[1]], [0, bias_est[2]], 'k')
ax.plot3D([0, sim.sensors.magnetometer.errormodel.bias_current[0]],
                           [0, sim.sensors.magnetometer.errormodel.bias_current[1]],
                           [0, sim.sensors.magnetometer.errormodel.bias_current[2]], 'g')

ax.legend(('True Bias', 'Estimated Bias', 'Data'))
ax.set_title('Magnetic Field [nT]')


with plt.rc_context(rc={'interactive': False}):
    plt.show()


