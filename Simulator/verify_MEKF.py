"""
Verify MEKF
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

import MEKF_cpp as MEKF
#-----------------------Clear Figures----------------------------
# clear figures
plt.close('all')

#-----------------Configuration / Parameters--------------------
tspan = np.array([0, 100])    # [sec]
L_cmd = np.zeros(3)			# initially command 0 torque


#----------------Initialize / Setup Workspace------------------
# setup sim
sim = Simulator(config)

# preallocate memory
T = np.arange(0, tspan[1]+config.tstep, config.tstep)
state_history = np.zeros((np.shape(T)[0], np.shape(sim.state)[0]))
state_history[0, :] = sim.state
B_body_history = np.zeros((np.shape(T)[0], 3))
command_history = np.zeros((np.shape(T)[0], 3))
DCM_history = np.zeros((np.shape(T)[0], 3, 3))
DCM_truth = np.zeros((np.shape(T)[0], 3, 3))
DCM_err = np.zeros(np.shape(T)[0])
#---------------------Propagate---------------------------
t = time.time()

# MEKF preallocation
xk = np.zeros(7,(np.shape(T)[0]))
Pk = np.zeros(6,6,(np.shape(T)[0])

W = 0.00000001 * np.identity(6) # needs to be updated
V = 0.0003 * np.identity(6) # needs to be updated

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
	B_ECI_pred = sim.environment.magfield_lookup(sim.state[0:3], sim.mag_order)
	
	V = np.array([B_ECI_pred/ np.linalg.norm(B_ECI_pred), S_ECI_pred])

	DCM_history[i, :, :] = triad.triad_ad(M.T, V.T)
	inter = conv.L(sim.state[3:7]) @ conv.R(sim.state[3:7]).T

	DCM_truth[i, :, :] = inter[1::, 1::]
	DCM_err[i] = np.linalg.norm(DCM_history[i, :, :] - DCM_truth[i, :, :].T)

    #--------------------MEKF-------------------------------
    
    xk[:,i+1] = MEKF.get_xk(xk[:,i],Pk[:,:,i],w_sensed,M,V,W,V,config.tstep)
    Pk[:,:,i+1] = MEKF.get_xk(xk[:,i],Pk[:,:,i],w_sensed,M,V,W,V,config.tstep)
    
    
    
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
plt.plot(T/3600, DCM_err)
plt.xlabel('time [hr]')
plt.title('DCM error')
plt.grid()

with plt.rc_context(rc={'interactive': False}):
	plt.show()

