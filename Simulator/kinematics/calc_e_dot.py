# -*- coding: utf-8 -*-

import numpy as np

def calc_e_dot(e,w):
    """
    Function: calc_e_dot

    Calculates derivative of euler angles.
    Input:
    e: euler angle vector (radians) np.array[3x1][float64]
    w: angular velocity vector (radians/s) np.array[3x1][float64]
    Output:
    e_dot: rate of change of eular angles np.array[1x3][float64]
    """
    #propagate actual dynamics
    tan_the = np.tan(e[1])
    if np.abs(tan_the) > 300: #make sure tan(theta) is well defined
        tan_the = 300*np.sign(tan_the)
    cos_the = np.cos(e[1])
    if np.abs(cos_the) < 10e-4: #make sure cos(theta) is well defined
        cos_the = 10e-4*np.sign(cos_the)
    A = np.array([[1, tan_the*np.sin(e[0]), tan_the*np.cos(e[0])],
                [0, np.cos(e[0]), -1*np.sin(e[0])],
                [0, np.sin(e[0])/cos_the, np.cos(e[0])/cos_the]])
    return np.dot(A,w)

#testing
e = np.transpose(np.array([np.pi/4,np.pi/4,np.pi/4]))
w = np.transpose(np.array([np.pi/4,np.pi/4,np.pi/4]))
print(calc_e_dot(e,w))

"""
CONVERT TO PYTHON

function e_dot   = calc_e_dot( e, w )
%INPUT: euler angles, angular velocity
%OUTPUT: rate of change of euler angles
%propagate actual dynamics
    tan_the = tan(e(2));
    if(abs(tan_the) > 300) %make sure tan(theta) is well defined
        tan_the = 300*sign(tan_the);
        disp('tan x')
    end
    cos_the = cos(e(2));
    if(abs(cos_the) < 10^-4) %make sure cos(theta) is well defined
        cos_the = 10^-4*sign(cos_the);
        disp('cos x')
    end
    A = [1, tan_the*sin(e(1)), tan_the*cos(e(1));
         0, cos(e(1)),              -1*sin(e(1));
         0, sin(e(1))/cos_the,  cos(e(1))/cos_the];
    e_dot = A*w;
end
"""
