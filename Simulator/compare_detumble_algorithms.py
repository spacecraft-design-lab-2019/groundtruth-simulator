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
max_dipoles = np.array([[8.8e-3], [1.373e-2], [8.2e-3]])


#----------------Initialize / Setup Workspace------------------
# setup sim
sim = Simulator(config)

# preallocate memory
T = np.arange(0, tspan[1]+config.tstep, config.tstep)
states_B_cross_bang_bang = np.zeros((np.shape(T)[0], np.shape(sim.state)[0]))
states_B_cross_directional = np.zeros((np.shape(T)[0], np.shape(sim.state)[0]))
states_B_dot = np.zeros((np.shape(T)[0], np.shape(sim.state)[0]))

B_ECI_history = np.zeros((np.shape(T)[0], 3))
B_body_history = np.zeros((np.shape(T)[0], 3))
command_history = np.zeros((np.shape(T)[0], 3))

#---------------------- B_dot --------------------------------------------------------------------------------------------
t = time.time()
sim = Simulator(config)

for i in range(len(T)):
	# Simulator
	sensors = sim.step(config.tstep, L_cmd)
	states_B_dot[i, :] = sim.state
	B_body = sim.debug_output[1]
	B_body_history[i,:] = B_body

	# command torque based on sensors (currently no noise addition 11/17)
	gain = 1 #.0143	#4e-2
	B_sensed = sensors[0:3]
	w_sensed = sensors[3:6]

	if i>0:
		B_1 = np.transpose(B_body_history[i-1,:])
		B_2 = np.transpose(B_body_history[i,:])
		B_dot = dcpp.get_B_dot(B_1,B_2,config.tstep)

		dipole = dcpp.detumble_B_dot_bang_bang(B_dot,max_dipoles)
		L_cmd = np.cross(np.squeeze(dipole), np.transpose(B_body)*1e-9)	 # Multiply by factor to put into SI units (Tesla)
		command_history[i,:] = np.transpose(dipole)

elapsed = time.time() - t
print(elapsed)

# #------------------------ B_cross_bang_bang -------------------------------------------------------------------------------
# t = time.time()
# sim2 = Simulator(config)
#
# for i in range(len(T)):
# 	# Simulator
# 	sensors = sim2.step(config.tstep, L_cmd)
# 	states_B_cross_bang_bang[i, :] = sim2.state
# 	B_body = sim2.debug_output[1]
#
# 	# command torque based on sensors (currently no noise addition 11/17)
# 	gain = 1 #.0143	#4e-2
# 	B_sensed = sensors[0:3]
# 	w_sensed = sensors[3:6]
# 	dipole = dcpp.detumble_B_cross_bang_bang(w_sensed, B_body, gain, max_dipoles)
# 	L_cmd = np.cross(np.squeeze(dipole), np.transpose(B_body)*1e-9)	 # Multiply by factor to put into SI units (Tesla)
# 	# print(i)
#
# elapsed = time.time() - t
# print(elapsed)
#
# #---------------------- B_cross_directional ----------------------------------------------------------------------------
# t = time.time()
# sim3 = Simulator(config)
#
# for i in range(len(T)):
# 	# Simulator
# 	sensors = sim3.step(config.tstep, L_cmd)
# 	states_B_cross_directional[i, :] = sim3.state
# 	B_body = sim3.debug_output[1]
#
# 	# command torque based on sensors (currently no noise addition 11/17)
# 	gain = 1 #.0143	#4e-2
# 	B_sensed = sensors[0:3]
# 	w_sensed = sensors[3:6]
# 	dipole = dcpp.detumble_B_cross_directional(w_sensed, B_sensed, gain, max_dipoles)
# 	L_cmd = np.cross(np.squeeze(dipole), np.transpose(B_body)*1e-9)	 # Multiply by factor to put into SI units (Tesla)
#
# elapsed = time.time() - t
# print(elapsed)
#
#


#------------------------Plot-----------------------------


plt.figure()
plt.plot(T/3600, np.linalg.norm(states_B_dot[:,10:13],axis=1))
plt.legend(('B_dot'))
plt.xlabel('time [hr]')
plt.title('angular velocity [rad/s]')
plt.grid()

plt.figure()
plt.plot(T/3600, states_B_dot[:,10])
plt.plot(T/3600, states_B_dot[:,11])
plt.plot(T/3600, states_B_dot[:,12])
plt.legend(('w_x', 'w_y', 'w_z'))
plt.xlabel('time [hr]')
plt.title('angular velocity [rad/s]')

plt.figure()
plt.plot(T/3600, states_B_dot[:,3])
plt.plot(T/3600, states_B_dot[:,4])
plt.plot(T/3600, states_B_dot[:,5])
plt.plot(T/3600, states_B_dot[:,6])
plt.legend(('q_1', 'q_2', 'q_3', 'q_4'))
plt.xlabel('time [hr]')
plt.title('Quaternion components')

plt.figure()
plt.plot(T/3600, B_body_history[:,0])
plt.plot(T/3600, B_body_history[:,1])
plt.plot(T/3600, B_body_history[:,2])
plt.legend(('B_x', 'B_y', 'B_z'))
plt.xlabel('time [hr]')
plt.title('Body magnetic field [uT]')

plt.figure()
plt.plot(T/3600, command_history[:,0])
plt.plot(T/3600, command_history[:,1])
plt.plot(T/3600, command_history[:,2])
plt.legend(('m_x', 'm_y', 'm_z'))
plt.xlabel('time [hr]')
plt.title('Commanded dipole')

with plt.rc_context(rc={'interactive': False}):
	plt.show()