# -*- coding: utf-8 -*-
import numpy as np
from kinematics import *
from dynamics import *
from initialize import *
from sciConsts import *

#-------------------------Simulator-------------------------------

def rk4_step(f, state, environment, h):
    """
    Function: rk4_step
        Uses rk4 to integrate a single step. (vectorized)

    Inputs:
        f:  derivative of state vector
        x: initial state object
        h:  step size (in seconds)
    Outputs:
        x1: updated state object
    """
    k1 = h * f(state, environment)

    x2 = copy(state).update_state(k1/2, h/2)
    k2 = h * f(x2, environment)

    x3 = copy(state).update_state(k2/2, h/2)
    k3 = h * f(x3, environment)

    x4 = copy(state).update_state(k3, h)
    k4 = h * f(x4, environment)

    xf = copy(state).update_state((k1 + 2*k2 + 2*k3 + k4)/6, h)







    x1 = copy(x).update_state(x.state()+k1/2)

    x.update_state(x.state()+k2/2)

    k3 = h * f(x.t+h/2,  x)
    x.update_state(x.state()+k3)

    k4 = h * f(x.t+h,    x)
    x1 = x.state() + (k1 + 2*k2 + 2*k3 + k4)/6

    t1 = x.t + h
    x.t = t1
    x.update_state(x1)
    return x


def calc_statedot(state, environment):
    """
    Function: state_dot
        Calculates the derivative of the state vector.

    Inputs:
        t: the current time
        s: the current state object
    Outputs:
        state_dot: the derivative of the state vector
    """

    #-----------------Calculate Environment --------------------------
    

    #----------------Calculate Accelerations/Torques------------------
    torque = np.zeros((3,))
    accel = np.zeros((3,))
    
    accel = accel + gravityPointMass(state.r, np.zeros((3,)), GM) # Earth Gravity
#    accel = accel + aeroDrag()
    
    torque = torque + gravityGradientTorque(state.r, state.I, GM)
    
    
    #-------------------Implement Control Law-------------------------
    
    

    #---------------------Kinematics----------------------------------
    q_dot = calc_q_dot(state.q, state.w)
    w_dot = calc_w_dot(state.w, torque, state.I)


    #---------------------Build Statedot------------------------------
    statedot = np.r_[state.v, q_dot, accel, w_dot]

    return statedot
