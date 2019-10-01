# -*- coding: utf-8 -*-
"""
Function: calc_q_dot

Calculates the derivative of the quaternion vector.
"""

def calc_q_dot(q, w):
    
    

"""
CONVERT TO PYTHON

function q_dot = calc_q_dot(q, w) 
%v2, try DAmico's notes
wx = w(1); wy = w(2); wz = w(3);
OMEGA = [   0,    wz, -1*wy, wx;
        -1*wz,     0,    wx, wy;
           wy, -1*wx,     0, wz;
        -1*wx, -1*wy, -1*wz,  0];
q_dot = 0.5*OMEGA*q;
end
"""