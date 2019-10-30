# -*- coding: utf-8 -*-
"""
Driver Script

Initializes and runs groundtruth simulator.

Pseudocode:
    -initialize with starting position & velocity

    LOOP:
        -calculate environment
        -calculate forces and torques (dynamics)
        -integrate (kinematics)
        -obtain new position/velocity

        -model sensors with noise

        -FEED to controller
        -take input from controller
        -add noise
        -account for control input in dynamics
"""

#------------ Initialize Workspace and Variables ------------

# Importing libraries
import numpy as np
import matplotlib.pyplot as plt
import julian
from datetime import datetime, timedelta
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv
from propagate_step import *
from constants import *


# Seed Start of Simulation with TLE
line1 = ('1 25544U 98067A   19302.69799672  .00001237  00000-0  29468-4 0  9996')
line2 = ('2 25544  51.6449  54.3108 0006382 201.3316 303.7490 15.50231045196128')


# Simulation Parameters
tstart = datetime(2019, 12, 30, 00, 00, 00)
tspan = np.array([0, 8640])    # [sec]
tstep = .1                     # [sec] - 10 Hz


# Setup SGP4
def sgp4_prog(S, t, microsecond = False):
    s = t.second + t.microsecond/1e6 if microsecond else t.second
    return S.propagate(t.year, t.month, t.day, t.hour, t.minute, s)

ISS_sgp4 = twoline2rv(line1, line2, wgs84)
r_i, v_i = sgp4_prog(ISS_sgp4, tstart)


# Initial Spacecraft State
q_i = np.array([1, 0, 0, 0])
w_i = np.array([.1, .5, -.3])
state_i = np.r_[r_i, q_i, v_i, w_i]


# Environment
mjd_start = julian.to_jd(tstart, 'mjd')
world = Environment(mjd_start=mjd_start)
struct = SpacecraftStructure()


# Initialize Variables
t = np.arange(0, tspan[1]+tstep, tstep)
state_history = np.zeros((np.shape(t)[0], np.shape(state_i)[0]))
state_history_sgp4 = np.zeros((np.shape(t)[0], 6))

state_history[0, :] = state_i
state_history_sgp4[0, :] = np.r_[r_i, v_i]


#------------------ Run Simulation -----------------------------

update = lambda t, state: calc_statedot(t, state, world, struct)
for idx in range(t.shape[0]-1):
    # Propagate
    state_history[idx+1, :] = rk4_step(update, t[idx], state_history[idx, :], tstep)

    # Normalize the Quaternion Vector
    state_history[idx+1, 3:7] = state_history[idx+1, 3:7]/np.linalg.norm(state_history[idx+1, 3:7])

    # SGP4
    dt = tstart + timedelta(seconds=t[idx])
    # for live feedback:
    print(dt)
    state_history_sgp4[idx+1, :] = np.array(sgp4_prog(ISS_sgp4, dt, microsecond = True)).flatten()


#-------------------- Plot------------------------------------

plt.figure()
plt.plot(state_history[:, 0], state_history[:, 1])
plt.xlabel("X_ECI")
plt.ylabel("Y_ECI")
plt.grid()
plt.show()

plt.figure()
plt.plot(t, state_history[:, 3:7])
plt.xlabel('time')
plt.ylabel('quaternions')
plt.grid()
plt.show()

plt.figure()
plt.plot(t, abs(state_history[:, 0:3] - state_history_sgp4[:, 0:3]))
plt.plot(t, abs(state_history[:, 7:10] - state_history_sgp4[:, 3:6]), hold='on')
plt.xlabel('time')
plt.ylabel('error')
plt.suptitle('Error against SGP4')
