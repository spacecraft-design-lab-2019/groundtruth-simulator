import numpy as np
from numpy import linalg as LA
def innovation(R,rN, rB, P,V,C):
    Q = np.zeros([6,6])
    Q[0:3,0:3] = R
    Q[3:6,3:6] = R

    
    z = rB - Q@rN

    S = C@P@C.T+V
    
    return z, S
