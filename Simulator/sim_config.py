# -*- coding: utf-8 -*-

import numpy as np
from datetime import datetime

#------------------- Configuration Parameters -------------------

# Seed Initial Position/Velocity with TLE - BEESAT-1
line1 = ('1 35933U 09051C   19315.45643387  .00000096  00000-0  32767-4 0  9991')
line2 = ('2 35933  98.6009 127.6424 0006914  92.0098 268.1890 14.56411486538102')

# Simulation Parameters
tstart = datetime(2019, 12, 30, 00, 00, 00)
tstep = .1                     # [sec] - 10 Hz

# Initial Spacecraft Attitude
q_i = np.array([1, 0, 0, 0])
w_i = np.array([.1, .5, -.3])


# Spacecraft Properties
I = np.array([[17,0,0],[0,18,0],[0,0,22]])
mass = 1.0 # kg

