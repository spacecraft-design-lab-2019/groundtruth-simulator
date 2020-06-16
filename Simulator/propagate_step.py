# -*- coding: utf-8 -*-
import numpy as np
import datetime, julian
from kinematics import *
from dynamics import *
from constants import Environment
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv

#-------------------------Simulator-------------------------------

def sgp4_step(line1, line2, mjd, model=wgs84):
    """
    Function: sgp4_step
        Returns position and velocity in ECI from SGP4 propagation.

    Inputs:
        line1:  first line of TLE (string)
        line2:  second line of TLE (string)
        mjd     time to propagate to (mjd)
        model:  Earth gravity model to use (default=wgs84)
    Outputs:
        rECI:   position in ECI
        vECI:   velocity in ECI
    """
    sgp4 = twoline2rv(line1, line2, model)
    dt = julian.from_jd(mjd, fmt='mjd')
    sec = dt.second + dt.microsecond/1e6
    return sgp4.propagate(dt.year, dt.month, dt.day, dt.hour, dt.minute, sec)


def rk4_step(f, t, state, h):
    """
    Function: rk4_step
        Uses rk4 to integrate a single step. (vectorized)

    Inputs:
        f:      function that differentiates state vector
        t:      current time
        state:  initial state object
        h:      step size
    Outputs:
        x1: updated state object
    """

    t1 = t + h/2.0
    t2 = t + h

    k1 = h * f(t,  state)
    k2 = h * f(t1, state + k1/2.0)
    k3 = h * f(t1, state + k2/2.0)
    k4 = h * f(t2, state + k3)

    x1 = state + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0

    return t2, x1


def calc_statedot(t, state, cmd, structure, environment, mag_order):
    """
    Function: state_dot
        Calculates the derivative of the state vector.

    Inputs:
        t:      current (simulation) time in seconds
        state:  current state object
        cmd:    commanded magnetic moment (Am^2 or Nm/T)
    Outputs:
        state_dot: the derivative of the state vector
    """

    r = state[0:3]
    q = state[3:7]
    v = state[7:10]
    w = state[10:13]
    T = state[13]


    #-----------------Calculate Environment --------------------------
    environment.update(t)

    
    #------------------Thermal Modelling------------------------------
    T_dot = structure.calc_Tdot(T, environment.isEclipse(r))

    
    #----------------Calculate Accelerations/Torques------------------
    torque = np.zeros(3)
    accel = np.zeros(3)

    adrag, mdrag = dragCalc(state, environment, structure)

    accel += gravityPointMass(r, environment.earth.GM)
    accel += gravityEarthJ2(r, environment.earth)
    accel += adrag

    torque += gravityGradientTorque(r, structure.I, environment.earth.GM)
    torque += mdrag


    #-------------------Implement Control Law-------------------------
    B_ECI = environment.magfield_lookup(r, mag_order)
    B_body = conv.quatrot(conv.conj(q),B_ECI)
    torque += np.cross(cmd, B_body);


    #---------------------Kinematics----------------------------------
    q_dot = calc_q_dot(q, w)
    w_dot = calc_w_dot(w, torque, structure.I)


    #---------------------Build Statedot------------------------------
    statedot = np.r_[v, q_dot, accel, w_dot, T_dot]

    return statedot
