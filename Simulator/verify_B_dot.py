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
B_ECI_history = np.zeros((np.shape(T)[0], 3))
B_body_history = np.zeros((np.shape(T)[0], 3))
command_history = np.zeros((np.shape(T)[0], 3))
#---------------------Propagate---------------------------
t = time.time()

for i, elapsed_t in enumerate(T[0:-1]):
	# Simulator
	sensors,B_ECI,B_body = sim.step(L_cmd, config.tstep)
	state_history[i+1, :] = sim.state
	B_ECI_history[i+1,:] = np.transpose(B_ECI)
	B_body_history[i+1,:] = np.transpose(B_body)
	command_history[i+1,:] = np.transpose(L_cmd)

	# command torque based on sensors (currently no noise addition 11/13)
	gain = .0143	#4e-2
	B_sensed = sensors[0:3]
	w_sensed = sensors[3:6]
	L_cmd = dcpp.detumble_B_cross(w_sensed, B_sensed, gain)

	# print(i)

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
plt.title('angular rates')
plt.grid()

plt.figure()
plt.plot(T/3600, B_ECI_history)
plt.xlabel('time [hr]')
plt.title('B_field_ECI [nT]')
plt.grid()

plt.figure()
plt.plot(T/3600, B_body_history)
plt.xlabel('time [hr]')
plt.title('B_field_body [nT]')
plt.grid()

plt.figure()
plt.plot(T/3600, command_history)
plt.xlabel('time [hr]')
plt.title('command history [Nm]')
plt.grid()

plt.figure()
plt.plot(T/3600, np.linalg.norm(state_history[:,10:13],axis=1))
plt.xlabel('time [hr]')
plt.title('angular velocity [rad/s]')
plt.grid()
with plt.rc_context(rc={'interactive': False}):
	plt.show()