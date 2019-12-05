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


''' currently taken from MEKF starter code - input real measurements/covariances here '''
rN1 = MEKF_inputs['rN1'] # unit vector measurement 2 in newtonian frame
rN2 = MEKF_inputs['rN2'] # unit vector measurement 2 in newtonian frame
rN = np.vstack((rN1,rN2))
whist = MEKF_inputs['whist'] # gyro measurement
rB2hist = MEKF_inputs['rB2hist']  # unit vector measurement 1 in body frame
rB1hist = MEKF_inputs['rB1hist'] # unit vector measurement 2 in body frame

dt = 0.1

# initial conditions:
rB1 = np.reshape(rB1hist[:,0],(3,1))
rB2 = np.reshape(rB2hist[:,0],(3,1))
Pk  = (10*math.pi/180)**2*np.eye(6) 
R0  = deterministic_ad.triad_ad(np.hstack((rB1,rB2)),np.hstack((rN1,rN2)))
q0  = DCM2q.DCM2quat(R0)  

# process and measurement noise (update with true values)
W = np.array(np.identity(6))*1e-9*0.3046
V = np.array(np.identity(6))
for ii in range(3):
    V(ii,ii) = 0.003
    V(ii+3,ii+3) = 0.0076

# initial conditions for MEKF
q0 = np.reshape(q0,(4,1))
Beta0 = np.reshape(np.array([0.1,0.1,0.1]),(3,1))
xk    = np.vstack((q0,Beta0))

# pre-allocate 
q_mekf  = np.zeros([4,1501])
b_mekf  = np.zeros([3,1501])

# initialize MEKF
q_mekf[0:4,0]  = np.reshape(q0,4)
b_mekf[0:3,0]  = np.reshape(Beta0,3)

# start sim
for ii in range(len(rB1hist.T)-1):
    # run predict step
    '''
    xn - predicted state (mu_k+1|k)
    Pn - predicted covariance (sigma_k+1|k)
    A  - linearized state transition matrix
    W  - Process noise covariance
    V  - Measurement noise covariance
    rN - Vector measurement in newtonian frame
    rB - Vector measurement in body frame
    '''
    # get measurement vector
    rB1 = np.reshape(rB1hist[:,ii+1],(3,1))
    rB2 = np.reshape(rB2hist[:,ii+1],(3,1))
    rB  = np.vstack((rB1,rB2))
    
    xn,A = predict(xk,whist[:,ii+1],dt)
    Pn   = A*Pk*np.transpose(A)+W
    
    # run measurement step
    y,R,C = measurement(xn[0:4],rN)
    
    # run innovation step to find z
    z,S = innovation(R, rN, rB, Pn, V, C)
    
    # Kalman Gain
    L = Pn @ np.transpose(C) @ LA.inv(S)
    
    # update step
    xk, Pk = update(L, z, xn, Pn, V, C)
    
    q_mekf[0:4,ii+1]  = np.reshape(xk[0:4],4)
    b_mekf[0:3,ii+1]  = np.reshape(xk[4:7],3)
