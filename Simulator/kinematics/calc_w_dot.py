# -*- coding: utf-8 -*-
import numpy as np

"""
Function: calc_w_dot

Uses Euler equations to calculate derivative of angular velocity vector.
"""

def calc_w_dot(w, torque, I):
    """
    Function: calc_w_dot
        Calculates the derivative of the angular velocity vector.
        
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
    S = np.array([0, -v[2], v[1]],
                 [v[2], 0, -v[0]],
                 [-v[1], v[0], 0])
    return S



#
#"""
#CONVERT TO PYTHON
#
#function w_dot = fcn(w, tor, I)
#
w_dot = -I\(Skew(w)*I*w-tor);
#
#end
#
#
#function S = Skew(v)
#S = [0 -v(3) v(2)
#    v(3) 0 -v(1)
#    -v(2) v(1) 0];
#end
#
#"""