# -*- coding: utf-8 -*-
"""
Initialization Script

Sets simulation parameters and initializes spacecraft. Run by driver-script
before beginning simulation.

NOTE: based on \Reference\Hridu Old Simulator\Model Scripts\initialize.m
"""
#Importing libraries
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

mu = 6.67408e-11*2e30

# Simulation Parameters
class SpacecraftParams():
    """
    A class to store all spacecraft parameters
    """
    def __init__(self,
                I = np.array([[17,0,0],[0,18,0],[0,0,22]]),
                r = np.array([0,1.5e11,0]),
                v = np.array([30000,0,0]),
                q = np.array([1,1,1,1])):
        self.I = I
        self.r = r
        self.v = v
        self.q = q

#testing SpacecraftParams
# s = SpacecraftParams()
# print(s.I)
# s.r = np.array([5,5,5])
# print(s.r)
# s.v[0] = 9
# print(s.v)

########## integrator #################
# I tried to plot a simple orbit (the moon around the earth)
# but it looks weird.
def calc_s_dot(t, s):
    con = -mu/(np.sqrt(s[0]**2 + s[1]**2 + s[2]**2)**3)
    return np.array([s[3],s[4],s[5],con*s[0],con*s[1],con*s[2]])

def one_step_integration(s, int_params, step):
    return int_params.integrate(step)

def driver(craft,t0,t1):
    s = np.concatenate((craft.r, craft.v), axis=0)
    t = np.linspace(t0, t1, 10000)
    y = np.zeros((len(t), len(s)))   # array for solution
    y[0, :] = s
    integrator = integrate.ode(calc_s_dot).set_integrator("dopri5")
    integrator.set_initial_value(s, t0)
    for i in range(1, t.size):
       y[i, :] = one_step_integration(y[i-1, :],integrator,t[i])
    plt.plot(y[:,1], y[:,0])
    plt.show()


driver(SpacecraftParams(),0,3.156e+7)


# Initial Spacecraft State Vector


# Mass Properties of Spacecraft


# Structure of Spacecraft


# Magnetorquers


# Sensors


# EKF Parameters


# Global Constants
