# -*- coding: utf-8 -*-
import numpy as np

#---------------------------Kinematics---------------------------------

#def calc_e_dot(e,w):
#    """
#    Function: calc_e_dot
#
#    Calculates derivative of euler angles.
#    Input:
#    e: euler angle vector (radians) np.array[3x1][float64]
#    w: angular velocity vector (radians/s) np.array[3x1][float64]
#    Output:
#    e_dot: rate of change of eular angles np.array[1x3][float64]
#    """
#    #propagate actual dynamics
#    tan_the = np.tan(e[1])
#    if np.abs(tan_the) > 300: #make sure tan(theta) is well defined
#        tan_the = 300*np.sign(tan_the)
#    cos_the = np.cos(e[1])
#    if np.abs(cos_the) < 10e-4: #make sure cos(theta) is well defined
#        cos_the = 10e-4*np.sign(cos_the)
#    A = np.array([[1, tan_the*np.sin(e[0]), tan_the*np.cos(e[0])],
#                [0, np.cos(e[0]), -1*np.sin(e[0])],
#                [0, np.sin(e[0])/cos_the, np.cos(e[0])/cos_the]])
#    return np.dot(A,w)


def calc_q_dot(q, w):
    """
    Function: calc_q_dot
        Calculates the derivative of the quaternion vector.
    
    Inputs:
        q: quaternion vector np.array[4x1][float64]
        w: angular velocity (radians/s) np.array[3x1][float64]
    Outputs:
        q_dot: rate of change of quaternion np.array[4x1][float64]
    """
    omega = np.array([[0, w[2], -1*w[1], w[0]],
                    [-1*w[2], 0, w[0], w[1]],
                    [w[1], -1*w[0], 0, w[2]],
                    [-1*w[0], -1*w[1], -1*w[2],0]])
    return 0.5*np.dot(omega,q)


def calc_w_dot(w, torque, I):
    """
    Function: calc_w_dot
        Calculates the derivative of the angular velocity vector using the
        Euler equations for rigid body dynamics.
        
    Inputs:
        w: angular velocity (radians/s) np.array[3x1][float64]
        torque: applied torque (Nm) np.array[3x1][float54]
    Ouputs:
        w_dot: derivative of angular velocity (rad/s/s) np.array[3x1][float64]
    """
    
    w_dot = np.linalg.solve(I, np.dot(skew(w), np.dot(I, w)) - torque)
    return w_dot
    

def skew(v):
    """
    Function: skew
        Calculates the skew matrix for cross-product calculation
        
    Inputs:
        
    """
    S = np.array([[0, -v[2], v[1]],
                 [v[2], 0, -v[0]],
                 [-v[1], v[0], 0]])
    return S


