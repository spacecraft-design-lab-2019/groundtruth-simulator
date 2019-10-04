# -*- coding: utf-8 -*-
import numpy as np
from kinematics import *
from dynamics import *
from initialize import *

import sys #temporary for debugging

#-------------------------Simulator-------------------------------

def rk4_step(f, x, h):
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
    test = SpacecraftState()
    k1 = h * f(x.t,      x).state()
    x.update_state(x.state()+k1/2)
    k2 = h * f(x.t+h/2,  x).state()
    x.update_state(x.state()+k2/2)
    k3 = h * f(x.t+h/2,  x).state()
    x.update_state(x.state()+k3)
    k4 = h * f(x.t+h,    x).state()
    x1 = x.state() + (k1 + 2*k2 + 2*k3 + k4)/6
    t1 = x.t + h
    x.t = t1
    x.update_state(x1)
    return x


def calc_statedot(t, state):
    """
    Function: state_dot
        Calculates the derivative of the state vector

    Inputs:
        t: the current time
        s: the current state object
    Outputs:
        state_dot: the derivative of the state vector
    """
    # fix: get I and other spacecraft params into the function some other way
    #I = np.array([[17,0,0],[0,18,0],[0,0,22]])
    GM = 3.986e5

    #pos = state[0:3]
    #q   = state[3:7]
    #vel = state[7:10]
    #w   = state[10:13]

    #-----------------Calculate Environment --------------------------


    #----------------Calculate Accelerations/Torques------------------
    torque = np.zeros((3,))
    accel = np.zeros((3,))
    #accel = accel + accelPointMass(pos, pos, GM)
    accel = accel + accelPointMass(state.r, state.r, GM)

    #---------------------Kinematics----------------------------------
    #q_dot = calc_q_dot(q, w)
    q_dot = calc_q_dot(state.q, state.w)
    #w_dot = calc_w_dot(w, torque, I)
    w_dot = calc_w_dot(state.w, torque, state.I)


    #---------------------Build Statedot------------------------------
    statedot = SpacecraftState(state.I,state.v,q_dot,accel,w_dot,state.t)

    #statedot = np.zeros(np.shape(state))
    #statedot[0:3]   = vel
    #statedot[3:7]   = q_dot
    #statedot[7:10]  = accel
    #statedot[10:13] = w_dot

    return statedot
