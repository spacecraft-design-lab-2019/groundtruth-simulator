# -*- coding: utf-8 -*-

"""
Verify detumble algorithms
"""

import pdb
import numpy as np
import matplotlib.pyplot as plt
import time
from sim_config import *
from simulation_step import simulation_step
from propagate_step import sgp4_step

# import GNC functions
import sys
sys.path.append('/home/eleboeuf/Documents/GNC')
import detumble_cpp as dcpp

#----------------Clear Figures---------------------------------
# clear figures
plt.close('all')

#----------------Initialize / Setup Workspace------------------
tspan = np.array([0, 600])    # [sec]
T = np.arange(0, tspan[1]+tstep, tstep)

# initialize classes
spacecraft = SpacecraftStructure(config.I, mass=config.mass)



#---------------------Initial State Vector---------------------
r_i, v_i = sgp4_step(line1, line2, tstart)
# pdb.set_trace()
state_i = np.r_[r_i, q_i, v_i, w_i]
state_history = np.zeros((np.shape(T)[0], np.shape(state_i)[0]))
state_history_sgp4 = np.zeros((np.shape(T)[0], 6))
state_history[0, :] = state_i
state_history_sgp4[0, :] = np.r_[r_i, v_i]
sim_state = {'state': state_i, 't': tstart}
L_cmd = np.zeros(3)			# initially command 0 torque


#---------------------Propagate---------------------------
t = time.time()
for i, elapsed_t in enumerate(T[0:-1]):
	# Simulator
	sensors, sim_state = simulation_step(L_cmd, sim_state)
	state_history[i+1, :] = sim_state['state']

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