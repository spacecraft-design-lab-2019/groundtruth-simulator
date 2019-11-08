# -*- coding: utf-8 -*-
import numpy as np
from kinematics import *
from dynamics import *
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv


#-------------------------Simulator-------------------------------

def sgp4_step(line1, line2, dt, model=wgs84):
    """
    Function: sgp4_step
        Returns position and velocity in ECI from SGP4 propagation.

    Inputs:
        line1:  first line of TLE (string)
        line2:  second line of TLE (string)
        dt:     time to propagate to (datetime)
        model:  Earth gravity model to use (default=wgs84)
    Outputs:
        rECI:   position in ECI
        vECI:   velocity in ECI
    """
    sgp4 = twoline2rv(line1, line2, model)
    sec = t.second + t.microsecond/1e6 if microsecond else t.second
    return sgp4.propagate(t.year, t.month, t.day, t.hour, t.minute, sec)


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


def calc_statedot(t, state, cmd, environment, structure):
    """
    Function: state_dot
        Calculates the derivative of the state vector.

    Inputs:
        t:      current time
        state:  current state object
        cmd:    controller input
    Outputs:
        state_dot: the derivative of the state vector
    """

    r = state[0:3]
    q = state[3:7]
    v = state[7:10]
    w = state[10:13]


    #-----------------Calculate Environment --------------------------
    mjd = environment.mjd_start + t/(24*60*60)


    #----------------Calculate Accelerations/Torques------------------
    torque = np.zeros(3)
    accel = np.zeros(3)

    adrag, mdrag = dragCalc(state, mjd, environment, structure)

    accel = accel + gravityPointMass(r, environment.earth.GM)
    accel = accel + gravityEarthJ2(r, environment.earth)
    accel = accel + adrag

    torque = torque + gravityGradientTorque(r, structure.I, environment.earth.GM)
    torque = torque + mdrag


    #-------------------Implement Control Law-------------------------
    # TO-DO: implement control
    torque = torque + cmd;


    #---------------------Kinematics----------------------------------
    q_dot = calc_q_dot(q, w)
    w_dot = calc_w_dot(w, torque, structure.I)


    #---------------------Build Statedot------------------------------
    statedot = np.r_[v, q_dot, accel, w_dot]

    return statedot
