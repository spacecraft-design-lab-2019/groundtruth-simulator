# -*- coding: utf-8 -*-

"""
Verify detumble algorithms
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import sim_config as config
from simulator import Simulator
import conversions as conv

# import GNC functions
import sys
sys.path.append('/home/eleboeuf/Documents/GNC')
import detumble_cpp as dcpp
import triad_cpp as triad
import sun_utils_cpp
import euler_cpp
#-----------------------Clear Figures----------------------------
# clear figures
plt.close('all')


#-----------------Configuration / Parameters--------------------
tspan = np.array([0, 1])    # [sec]
L_cmd = np.zeros(3)			# initially command 0 torque


#----------------Initialize / Setup Workspace------------------
# setup sim
sim = Simulator(config)
max_dipoles = np.array([[8.8e-3], [1.373e-2], [8.2e-3]]) # TODO: Add dipole saturation to spacecraft properties

# preallocate memory
T = np.arange(0, tspan[1]+config.tstep, config.tstep)
state_history = np.zeros((np.shape(T)[0], np.shape(sim.state)[0]))
state_history[0, :] = sim.state
B_body_history = np.zeros((np.shape(T)[0], 3))
command_history = np.zeros((np.shape(T)[0], 3))
DCM_history = np.zeros((np.shape(T)[0], 3, 3))
DCM_truth = np.zeros((np.shape(T)[0], 3, 3))
#---------------------Propagate---------------------------
t = time.time()

for i, elapsed_t in enumerate(T[0:-1]):
	# Simulator
	sensors = sim.step(config.tstep, L_cmd)
	state_history[i+1, :] = sim.state
	
	command_history[i+1,:] = np.transpose(L_cmd)

	# command torque based on sensors (currently no noise addition 11/17)
	B_sensed = sensors[0:3]
	w_sensed = sensors[3:6]
	S_sensed = sensors[6:9]
	B_body_history[i+1,:] = np.transpose(B_sensed)
	M = np.array([B_sensed / np.linalg.norm(B_sensed), S_sensed])

	S_ECI_pred = sun_utils_cpp.sat_sun_vect(sim.state[0:3], sim.MJD) 
	S_ECI_pred = S_ECI_pred / np.linalg.norm(S_ECI_pred)
	B_ECI_pred = sim.environment.magfield_lookup(sim.state[0:3])
	
	V = np.array([B_ECI_pred/ np.linalg.norm(B_ECI_pred), S_ECI_pred])

	DCM_history[i, :, :] = triad.triad_ad(M.T, V.T)
	inter = conv.L(sim.state[3:7]) @ conv.R(sim.state[3:7]).T
	inter2 = euler_cpp.Lq(sim.state[3:7]) @ euler_cpp.Rq(sim.state[3:7]).T
	if i == 0:
		print(inter[1::, 1::])
		print(inter2)
		print(triad.triad_ad(M.T, V.T))
	DCM_truth[i, :, :] = inter[1::, 1::]

	#--------------------B_cross---------------------------
	gain_B_cross = .0143  # 4e-2
	L_cmd = dcpp.detumble_B_cross(w_sensed, B_sensed, gain_B_cross)

	#--------------------B_dot-----------------------------
	# if i>1:
	# 	gain_B_dot = 5e-6
	# 	B_dot = dcpp.get_B_dot(np.transpose(B_body_history[i,:]), B_sensed, config.tstep)
	# 	dipole = dcpp.detumble_B_dot_bang_bang(B_dot, max_dipoles)
	# 	bang_bang_gain = 1e-9  # 5e-6
	# 	L_cmd = np.cross(dipole,B_sensed)


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

plt.figure()
plt.plot(T/3600, DCM_history[:, 0, 1] - DCM_truth[:, 0, 1])
# plt.plot(T/3600, )

with plt.rc_context(rc={'interactive': False}):
	plt.show()

