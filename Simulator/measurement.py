import sys
import numpy as np

def SkewSymmetric(a):
    a1 = a[0]
    a2 = a[1]
    a3 = a[2]
    a_skew = np.array([[0,-a3,a2],[a3,0,-a1],[-a2,a1,0]])
    return a_skew
def quat2DCM(quat):
    '''
    Takes in a quaternion converts to a direction cosine matrix
    Inputs:
        q - quaternion, scalar first, (4x1)
        !!!!!!!! Why is this scalar first, but q_dot function uses scalar last!!!!!!-Paul
    Outputs:
        DCM - 3x3 direction cosine matrix
    TODO:
        - clarify quaternion convention
        - clarify quaternion rotation direction convention (eci2attitude, attitude2eci?)
        - clarify DCM rotation direction convention
    '''
    quat = np.reshape(quat,4)
    q1 = quat[0]
    q2 = quat[1]
    q3 = quat[2]
    q4 = quat[3]
    q_vec = np.array([[q2],[q3],[q4]])
    # calculate skew symmetric matrix Q
    Q = np.array([[0,-q4,q3],[q4,0,-q2],[-q3,q2,0]])
    # calculate DCM
    DCM = (q1**2-np.transpose(q_vec)@q_vec)*np.identity(3)-2*q1*Q+2*q_vec@np.transpose(q_vec)

    return DCM


def measurement(q,rN):
    ''' 
    R = Rotation from N frame to B frame
    '''

    R = quat2DCM(q)

    rB1 = (R@rN[0:3])
    rB2 = (R@rN[3:6])
    y = np.append(rB1,rB2)
    C = np.zeros([6,6])
    C[0:3,0:3] = 2*SkewSymmetric(rB1)
    C[3:6,0:3] = 2*SkewSymmetric(rB2)
    
    return y,R,C
