# -*- coding: utf-8 -*-
"""
Initialization Script

Sets simulation parameters and initializes spacecraft. Run by driver-script
before beginning simulation.

"""

#-------------------------Setup---------------------------------

# Importing libraries
import numpy as np
import julian
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv


# Example ISS (Zarya) TLE
line1 = ('1 25544U 98067A   19302.69799672  .00001237  00000-0  29468-4 0  9996')
line2 = ('2 25544  51.6449  54.3108 0006382 201.3316 303.7490 15.50231045196128')


# Simulation Parameters
tstart = datetime(2019, 10, 30, 00, 00, 00)
tspan = np.array([0, 8640])    # [sec]
tstep = .1                     # [sec] - 10 Hz


# Initialize
ISS_sgp4 = twoline2rv(line1, line2, wgs84)
r_i, v_i = ISS_sgp4.propagate(tstart.year, tstart.month, tstart.day, tstart.hour, tstart.minute, tstart.second)
q_i = np.array([1, 0, 0, 0])
w_i = np.array([.1, .5, -.3])
state_i = np.r_[r_i, q_i, v_i, w_i]

mjd_start = julian.to_jd(tstart, 'mjd')


#--------------------------------------------------------------