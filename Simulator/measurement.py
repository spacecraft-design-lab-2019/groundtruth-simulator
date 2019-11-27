import sys
import numpy as np

sys.path.append('../../euler/') 
import Euler

def SkewSymmetric(a):
    a1 = a[0]
    a2 = a[1]
    a3 = a[2]
    a_skew = np.array([[0,-a3,a2],[a3,0,-a1],[-a2,a1,0]])
    return a_skew

def measurement(q,rN):
    ''' 
    R = Rotation from N frame to B frame
    '''

    R = Euler.quat2DCM(q)

    rB1 = (R@rN[0:3])
    rB2 = (R@rN[3:6])
    y = np.append(rB1,rB2)
    C = np.zeros([6,6])
    C[0:3,0:3] = 2*SkewSymmetric(rB1)
    C[3:6,0:3] = 2*SkewSymmetric(rB2)
    
    return y,R,C