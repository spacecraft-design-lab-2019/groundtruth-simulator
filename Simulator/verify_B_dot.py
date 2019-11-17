# -*- coding: utf-8 -*-

"""
Verify detumble algorithms
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import sim_config as config
from simulator import Simulator


# import GNC functions
import sys
sys.path.append('/home/eleboeuf/Documents/GNC')
import detumble_cpp as dcpp

#-----------------------Clear Figures----------------------------
# clear figures
plt.close('all')


#-----------------Configuration / Parameters--------------------
tspan = np.array([0, 600])    # [sec]
L_cmd = np.zeros(3)			# initially command 0 torque


#----------------Initialize / Setup Workspace------------------
# setup sim
sim = Simulator(config)

# preallocate memory
T = np.arange(0, tspan[1]+config.tstep, config.tstep)
state_history = np.zeros((np.shape(T)[0], np.shape(sim.state)[0]))
state_history[0, :] = sim.state

#---------------------Propagate---------------------------
t = time.time()
for i, elapsed_t in enumerate(T[0:-1]):
	# Simulator
	sensors = sim.step(L_cmd, config.tstep)
	state_history[i+1, :] = sim.state

	# command torque based on sensors (check simulation_step.py for noise addition, currently no noise addition 11/13)
	gain = 4e-12	#4e-2
	B_sensed = sensors[0:3]
	w_sensed = sensors[3:6]
	L_cmd = dcpp.detumble_B_cross(w_sensed, B_sensed, gain)
	print(i)

elapsed = time.time() - t
print(elapsed)

#------------------------Plot-----------------------------
plt.figure()
plt.plot(T/3600, state_history[:,3:7])
plt.xlabel('time [hr]')
plt.ylabel('quaternions')
plt.grid()

plt.figure()
plt.plot(T/3600, state_history[:,10:13])
plt.xlabel('time [hr]')
plt.ylabel('quaternions')
plt.grid()

plt.figure()
plt.plot(T/3600, np.linalg.norm(state_history[:,10:13],axis=1))
plt.xlabel('time [hr]')
plt.ylabel('angular velocity [rad/s]')
plt.grid()
with plt.rc_context(rc={'interactive': False}):
	plt.show()