# -*- coding: utf-8 -*-

import numpy as np
from datetime import datetime

#------------------- Configuration Parameters -------------------

# Seed Initial Position/Velocity with TLE
line1 = ('1 25544U 98067A   19302.69799672  .00001237  00000-0  29468-4 0  9996')
line2 = ('2 25544  51.6449  54.3108 0006382 201.3316 303.7490 15.50231045196128')

# Simulation Parameters
tstart = datetime(2019, 12, 30, 00, 00, 00)
tspan = np.array([0, 8640])    # [sec]
tstep = .1                     # [sec] - 10 Hz

# Initial Spacecraft Attitude
q_i = np.array([1, 0, 0, 0])
w_i = np.array([.1, .5, -.3])


# Spacecraft Properties
I = np.array([[17,0,0],[0,18,0],[0,0,22]]) 

