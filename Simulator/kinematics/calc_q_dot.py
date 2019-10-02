# -*- coding: utf-8 -*-
import numpy as np

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

##testing
#q = np.transpose(np.array([np.pi/4,np.pi/4,np.pi/4,np.pi/4]))
#w = np.transpose(np.array([np.pi/4,np.pi/4,np.pi/4]))
#print(calc_q_dot(q,w))
#
#"""
#CONVERT TO PYTHON
#
#function q_dot = calc_q_dot(q, w)
#%v2, try DAmico's notes
#wx = w(1); wy = w(2); wz = w(3);
#OMEGA = [   0,    wz, -1*wy, wx;
#        -1*wz,     0,    wx, wy;
#           wy, -1*wx,     0, wz;
#        -1*wx, -1*wy, -1*wz,  0];
#q_dot = 0.5*OMEGA*q;
#end
#"""
