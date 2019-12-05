import numpy as np
import math
from numpy import linalg as LA
# import from other MEKF files
from predict import predict
from measurement import measurement
from innovation import innovation
from update import update
# import functions from other filepaths
import sys
sys.path.append('/home/eleboeuf/Documents/GNC')

# for sim
import time
import sim_config as config
from simulator import Simulator
import conversions as conv

# import GNC functions

import detumble_cpp as dcpp
import triad_cpp as triad
import sun_utils_cpp
import euler_cpp


dt = 0.1

# initial conditions:
Pk  = (10*math.pi/180)**2*np.eye(6) 


# process and measurement noise (update with true values)
W = np.array(np.eye(6))*1e-9*0.3046
V = np.array(np.eye(6))
for ii in range(3):
    V(ii,ii) = 0.003
    V(ii+3,ii+3) = 0.0076

# initial conditions for MEKF
Beta0 = np.reshape(np.array([0.1,0.1,0.1]),(3,1))

# pre-allocate 
q_mekf  = np.zeros([4,1501])
b_mekf  = np.zeros([3,1501])

# initialize MEKF
q_mekf[0:4,0]  = np.reshape(q0,4)
b_mekf[0:3,0]  = np.reshape(Beta0,3)

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
	rB = np.array([B_sensed / np.linalg.norm(B_sensed), S_sensed])

	S_ECI_pred = sun_utils_cpp.sat_sun_vect(sim.state[0:3], sim.MJD) 
	S_ECI_pred = S_ECI_pred / np.linalg.norm(S_ECI_pred)
	B_ECI_pred = sim.environment.magfield_lookup(sim.state[0:3])
	
	rN = np.array([B_ECI_pred/ np.linalg.norm(B_ECI_pred), S_ECI_pred])

	DCM_history[i, :, :] = triad.triad_ad(rB.T, rN.T)
	inter = conv.L(sim.state[3:7]) @ conv.R(sim.state[3:7]).T
	inter2 = euler_cpp.Lq(sim.state[3:7]) @ euler_cpp.Rq(sim.state[3:7]).T
    
    '''
    xn - predicted state (mu_k+1|k)
    Pn - predicted covariance (sigma_k+1|k)
    A  - linearized state transition matrix
    W  - Process noise covariance
    V  - Measurement noise covariance
    rN - Vector measurement in newtonian frame
    rB - Vector measurement in body frame
    '''
    
    #rB = np.vstack((M[0:3],M[3:6]))
    #rN = np.vstack((V[0:3],V[3:6]))
    q0 = DCM2q.DCM2q(DCM_history[i,:,:])
    xk    = np.vstack((q0,Beta0))

    
    
    xn,A = predict(xk,w_sensed,dt)
    Pn   = A*Pk*np.transpose(A)+W
    
    # run measurement step
    y,R,C = measurement(xn[0:4],rN)
    
    # run innovation step to find z
    z,S = innovation(R, rN, rB, Pn, V, C)
    
    # Kalman Gain
    L = Pn @ np.transpose(C) @ LA.inv(S)
    
    # update step
    xk, Pk = update(L, z, xn, Pn, V, C)
    
    q_mekf[0:4,i+1]  = np.reshape(xk[0:4],4)
    b_mekf[0:3,i+1]  = np.reshape(xk[4:7],3)
