# -*- coding: utf-8 -*-
import numpy as np
from kinematics import *
from dynamics import *
# from constants import *

#-------------------------Simulator-------------------------------

def rk4_step(f, t, state, h):
    """
    Function: rk4_step
        Uses rk4 to integrate a single step. (vectorized)

    Inputs:
        f:      function that differentiates state vector
        t:      current time
        state:  initial state object
        h:      step size (in seconds)
    Outputs:
        x1: updated state object
    """

    k1 = h * f(t,       state)
    k2 = h * f(t + h/2, state + k1/2)
    k3 = h * f(t + h/2, state + k2/2)
    k4 = h * f(t + h,   state)

    return state + (k1 + 2*k2 + 2*k3 + k4)/6


def calc_statedot(t, state, environment, structure):
    """
    Function: state_dot
        Calculates the derivative of the state vector.

    Inputs:
        t:      the current time
        state:  the current state object
    Outputs:
        state_dot: the derivative of the state vector
    """

    r = state[0:3]
    q = state[3:7]
    v = state[7:10]
    w = state[10:13]


    #-----------------Calculate Environment --------------------------
    

    #----------------Calculate Accelerations/Torques------------------
    torque = np.zeros((3,))
    accel = np.zeros((3,))
    
    accel = accel + gravityPointMass(r, np.zeros((3,)), environment.earth.GM)
    accel = accel + gravityEarthJ2(r, environment.earth.GM, environment.earth.J2, environment.earth.radius)

    torque = torque + gravityGradientTorque(r, structure.I, environment.earth.GM)
    
    
    #-------------------Implement Control Law-------------------------
    
    

    #---------------------Kinematics----------------------------------
    q_dot = calc_q_dot(q, w)
    w_dot = calc_w_dot(w, torque, structure.I)


    #---------------------Build Statedot------------------------------
    statedot = np.r_[v, q_dot, accel, w_dot]

    return statedot
