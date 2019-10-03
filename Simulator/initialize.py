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

mu = 3.986e5

# Simulation Parameters
class SpacecraftParams():
    """
    A class to store all spacecraft parameters
    """
    def __init__(self,
                I = np.array([[17,0,0],[0,18,0],[0,0,22]]),
                r = np.array([363104000,0,0]),
                v = np.array([0,0.033132,0]),
                q = np.array([1,1,1,1])):
        self.I = I
        self.r = r
        self.v = v
        self.q = q

#testing SpacecraftParams
s = SpacecraftParams()
print(s.I)
s.r = np.array([5,5,5])
print(s.r)
s.v[0] = 9
print(s.v)

########## integrator #################
def rk4_step(f, t0, x0, h):
    """
    Function: rk4_step
        Uses rk4 to integrate a single step
        
    Inputs:
        f:  function we are integrating
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

# I tried to plot a simple orbit (the moon around the earth)
# but it looks weird.
def calc_s_dot(t, s):
    con = -mu/(np.sqrt(s[0]**2 + s[1]**2 + s[2]**2)**3)
    return np.array([s[3],s[4],s[5],con*s[0],con*s[1],con*s[2]])

def one_step_integration(s, int_params, step):
    return int_params.integrate(step)

def driver(craft,t0,t1):
    s = np.concatenate((craft.r, craft.v), axis=0)
    t = np.linspace(t0, t1, 100)
    y = np.zeros((len(t), len(s)))   # array for solution
    y[0, :] = s
    integrator = integrate.ode(calc_s_dot).set_integrator("dopri5")
    integrator.set_initial_value(s, t0)
    for i in range(1, t.size):
       y[i, :] = one_step_integration(y[i-1, :],integrator,t[i])
    plt.plot(y[:,1], y[:,0])
    plt.show()


driver(SpacecraftParams(),1,2332800)


# Initial Spacecraft State Vector


# Mass Properties of Spacecraft


# Structure of Spacecraft


# Magnetorquers


# Sensors


# EKF Parameters


# Global Constants
