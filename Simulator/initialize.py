# -*- coding: utf-8 -*-
"""
Initialization Script

Sets simulation parameters and initializes spacecraft. Run by driver-script
before beginning simulation.

NOTE: based on \Reference\Hridu Old Simulator\Model Scripts\initialize.m
"""
#Importing libraries
import numpy as np
from scipy.integrate import quad

GM = 3.986e5

# Simulation Parameters
class SpacecraftParams():
    """
    A class to store all spacecraft parameters
    """
    def __init__(self,
                I = np.array([[17,0,0],[0,18,0],[0,0,22]]),
                r = np.array([1,1,1]),
                v = np.array([1,1,1]),
                q = np.array([1,1,1,1])):
        self.I = I
        self.r = r
        self.v = v
        self.q = q

def driver(craft):
    s = np.concatenate((craft.r, craft.v), axis=0)
    s_dot = calc_s_dot(s,GM)
    quad()

def one_step_integrator(s,func,delta_t):
    t0 = 0
    y = np.zeros((1, len(s)))   # array for solution
    r = integrate.ode(func).set_integrator("dopri5")
    r.set_initial_value(s, t0)
    y[0, :] = r.integrate(delta_t) # get one more value, add it to the array
    if not r.successful():
        raise RuntimeError("Could not integrate")
    return y

def calc_s_dot(t,s):
    con = -mu/(np.sqrt(s[0]**2 + s[1]**2 + s[2]**2)**3)
    return np.array([s[3],s[4],s[5],con*s[0],con*s[1],con*s[2]])
    


def orbit_integrator:


#testing SpacecraftParams
s = SpacecraftParams()
print(s.I)
s.r = np.array([5,5,5])
print(s.r)
s.v[0] = 9
print(s.v)
# Initial Spacecraft State Vector


# Mass Properties of Spacecraft


# Structure of Spacecraft


# Magnetorquers


# Sensors


# EKF Parameters


# Global Constants
