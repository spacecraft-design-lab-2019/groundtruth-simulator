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
#this works but we need to change from a dummy function
def vdp1(t, y):
    return np.array([y[1], (1 - y[0]**2)*y[1] - y[0]])

def integrat(s, inte, step):
    return inte.integrate(step)
t = np.linspace(0,2332800,2332800/10000)

t0, t1 = 0, 20                # start and end
t = np.linspace(t0, t1, 1000)  # the points of evaluation of solution
y0 = [2, 0]                   # initial value
y = np.zeros((len(t), len(y0)))   # array for solution
y[0, :] = y0
integrator = integrate.ode(vdp1).set_integrator("dopri5")
integrator.set_initial_value(y0, t0)
for i in range(1, t.size):
   y[i, :] = integrat(y[i-1, :],integrator,t[i]) # get one more value, add it to the array

plt.plot(t, y)
plt.show()



#this doesn't work--still debugging
def driver(craft):
    s = np.concatenate((craft.r, craft.v), axis=0)
    s_dot = calc_s_dot(s,GM)
    quad()

def one_step_integrator(s,func,delta_t):
    t0 = 0
    r = integrate.ode(func).set_integrator("dopri5")
    r.set_initial_value(s, t0)
    y = r.integrate(delta_t)
    if not r.successful():
        raise RuntimeError("Could not integrate")
    return y

def calc_s_dot(t,s):
    con = -mu/(np.sqrt(s[0]**2 + s[1]**2 + s[2]**2)**3)
    return np.array([s[3],s[4],s[5],con*s[0],con*s[1],con*s[2]])


t = np.linspace(0,2332800,2332800/10000)
res = np.zeros((len(t), len(s)))
res[0, :] = s
for i in range(1, t.size):
    s = one_step_integrator(s,calc_s_dot,2332800/10000)
    res[i, :] = s

plt.plot(res[0],res[1])

def orbit_integrator:



# Initial Spacecraft State Vector


# Mass Properties of Spacecraft


# Structure of Spacecraft


# Magnetorquers


# Sensors


# EKF Parameters


# Global Constants
