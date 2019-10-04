# -*- coding: utf-8 -*-
import numpy as np
from kinematics import *
from dynamics import *

#-------------------------Simulator-------------------------------

def rk4_step(f, t0, x0, h):
    """
    Function: rk4_step
        Uses rk4 to integrate a single step. (vectorized)
        
    Inputs:
        f:  derivative of state vector
        t0: initial time
        x0: initial state vector
        h:  step size (in seconds)
    Outputs:
        t1: updated time
        x1: updated state vector
    """
    
    k1 = h * f(t0,      x0)
    k2 = h * f(t0+h/2,  x0+k1/2)
    k3 = h * f(t0+h/2,  x0+k2/2)
    k4 = h * f(t0+h,    x0+k3)
    x1 = x0 + (k1 + 2*k2 + 2*k3 + k4)/6
    t1 = t0 + h
    return t1, x1


def calc_statedot(t, state):
    """
    Function: state_dot
        Calculates the derivative of the state vector
        
    Inputs:
        t: the current time
        s: the current state vector
    Outputs:
        state_dot: the derivative of the state vector
    """
    
    # fix: get I and other spacecraft params into the function some other way
    I = np.array([[17,0,0],[0,18,0],[0,0,22]])
    GM = 3.986e5

    pos = state[0:3]
    q   = state[3:7]
    vel = state[7:10]
    w   = state[10:13]
    
    #-----------------Calculate Environment --------------------------
    
    
    #----------------Calculate Accelerations/Torques------------------
    torque = np.zeros((3,))
    accel = np.zeros((3,))
    accel = accel + accelPointMass(pos, pos, GM)
    
    
    #---------------------Kinematics----------------------------------
    q_dot = calc_q_dot(q, w)
    w_dot = calc_w_dot(w, torque, I)
    
    
    #---------------------Build Statedot------------------------------
    statedot = np.zeros(np.shape(state))
    statedot[0:3]   = vel
    statedot[3:7]   = q_dot
    statedot[7:10]  = accel
    statedot[10:13] = w_dot
    
    return statedot