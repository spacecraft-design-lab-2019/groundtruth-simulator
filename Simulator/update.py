import numpy as np
import math
from numpy import linalg as LA
def quatMult(q1,q2):
    qn = np.zeros(4)
    qn[0] = q1[0]*q2[0]-np.matmul(q1[1:4].T,q2[1:4])
    
    qn_temp = q1[0]*q2[1:4]+q2[0]*q1[1:4]+np.reshape(np.cross(q1[1:4].T,q2[1:4].T),(3,1))
    
    qn[1:4] = np.reshape(qn_temp,3)
    qn = qn / LA.norm(qn)
    # propagate dynamics forward
    return np.reshape(qn,(4,1))

def update(L, z, xn, Pn, V, C):
    
    dx = L@z
    dphi = dx[0:3]
    dbeta = dx[3:6]

    temp = math.sqrt(1-np.reshape(dphi.T@dphi,1))
    temp_q = np.vstack((temp,dphi))
    xn = np.reshape(xn,(7,1))
    q = quatMult(xn[0:4],temp_q)

    b = xn[4:7] + dbeta
    b = np.reshape(b,(3,1))
    
    x = np.vstack((q,b))
    P = (np.identity(6)-L*C)*Pn*np.transpose(np.identity(6)-L*C)+L*V*np.transpose(L)
    
    return x, P