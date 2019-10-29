# -*- coding: utf-8 -*-
import numpy as np
import pdb
import julian
from datetime import datetime
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv


#-------------------------Setup---------------------------------

# Example ISS (Zarya) TLE
line1 = ('1 25544U 98067A   19302.69799672  .00001237  00000-0  29468-4 0  9996')
line2 = ('2 25544  51.6449  54.3108 0006382 201.3316 303.7490 15.50231045196128')

# Start Time
t_i = datetime(2019, 10, 30, 00, 00, 00)

# Initial Position/Velocity
ISS = twoline2rv(line1, line2, wgs84)
r_i, v_i = ISS.propagate(t_i.year, t_i.month, t_i.day, t_i.hour, t_i.minute, t_i.second)

from inititalize import *


#------------------ Run Simulation -----------------------------

world = Environment(mjd_start=mjd_start)
struct = SpacecraftStructure()

t = np.arange(0, tspan[1]+tstep, tstep)
state_history = np.zeros((np.shape(t)[0], np.shape(state_i)[0]))
state_history[0, :] = state_i

for idx in range(t.shape[0]-1):
    # Integrate
    state_history[idx+1, :] = rk4_step(lambda t, state: calc_statedot(t, state, world, struct), t[idx], state_history[idx, :], tstep)

    # Normalize the Quaternion Vector
    state_history[idx+1, 3:7] = state_history[idx+1, 3:7]/np.linalg.norm(state_history[idx+1, 3:7])