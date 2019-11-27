import numpy as np
import math
from numpy import linalg as LA
# import from other files
from predict import predict
from measurement import measurement
from innovation import innovation
from update import update
# import functions from other filepaths
import sys
sys.path.append('../../TRIAD/py_funcs/')
import deterministic_ad
sys.path.append('../../util_funcs/py_funcs/')
import DCM2q 
# import data from sample mat file
from scipy.io import loadmat
#from scipy.io import savemat

MEKF_inputs = loadmat('mekf_inputs.mat')
#MEKF_truth  = loadmat('mekf_truth.mat')
#b_true = MEKF_truth['btrue']
#q_true = MEKF_truth['qtrue']
rN1 = MEKF_inputs['rN1'] # need updated
rN2 = MEKF_inputs['rN2'] # need updated
rN = np.vstack((rN1,rN2))
whist = MEKF_inputs['whist'] # from state machine
rB2hist = MEKF_inputs['rB2hist'] # from state machine???
rB1hist = MEKF_inputs['rB1hist']
W = MEKF_inputs['W'] # need updated
V = MEKF_inputs['V'] # need updated
dt = MEKF_inputs['dt']

# initial conditions:
rB1 = np.reshape(rB1hist[:,0],(3,1))
rB2 = np.reshape(rB2hist[:,0],(3,1))
Pk  = (10*math.pi/180)**2*np.eye(6) 
R0  = deterministic_ad.triad_ad(np.hstack((rB1,rB2)),np.hstack((rN1,rN2)))
q0  = DCM2q.DCM2quat(R0)  

#q0 = np.array([0.2737,0.0438,-0.9485,0.1535])
q0 = np.reshape(q0,(4,1))

# initial conditions for MEKF
Beta0 = np.reshape(np.array([0.1,0.1,0.1]),(3,1))
xk    = np.vstack((q0,Beta0))

# pre-allocate 
q_mekf  = np.zeros([4,1501])
b_mekf  = np.zeros([3,1501])
q_model = np.zeros([4,1501])

# initialize
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
# save to mat file to process in matlab
#savemat('output_mekf.mat', {
#        'q_mekf': q_mekf,
#        'b_mekf': b_mekf,
#        'q_model': q_model
#    })