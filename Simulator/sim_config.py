# -*- coding: utf-8 -*-

import numpy as np
from datetime import datetime
import os
import sys
import julian
from math import sqrt

#------------------- Configuration Parameters -------------------

# Seed Initial Position/Velocity with TLE - BEESAT-1
# (optional) - can instead replace this with r_i, v_i as np.array(3)
line1 = ('1 35933U 09051C   19315.45643387  .00000096  00000-0  32767-4 0  9991')
line2 = ('2 35933  98.6009 127.6424 0006914  92.0098 268.1890 14.56411486538102')

# Simulation Parameters
tstart = datetime(2019, 12, 30, 00, 00, 00)
tstep = .1                     # [sec] - 1 Hz
MJDstart = julian.to_jd(tstart, fmt='mjd')

# Initial Spacecraft State
q_i = np.array([sqrt(4.0)/4.0, sqrt(4.0)/4.0, sqrt(4.0)/4.0, sqrt(4.0)/4.0])    # quaternion
w_i = np.array([.03, .03, .03])   # radians/sec
T_i = 283 # Kelvin

# Spacecraft Properties
I = np.array([[1.959e-4, 2016.333e-9, 269.176e-9],[2016.333e-9, 1.999e-4, 2318.659e-9],[269.176e-9, 2318.659e-9, 1.064e-4]]) # kg*m^2
mass = 177.808e-3 # kg
thermal_properties = {
	"heat_capacitance" : 1000, # J/kgK
	"absorptivity" : 0.95,
	"emmisivity" : 0.89,
	"albedo" : 0.37
}

# Sensor Parameters
gyro_params = {
	"scalefactor" : 0.002,
	"crossaxis_sensitivity" : 0.02,
	"b" : sqrt(.0022),  # initialize with np.random and then don't change ever
	"cov" : 0.000000000694444,
	"random_walk_cov" : 0.00 	# will be mulitplied to np.random.rand(3)
}

mag_params = {
	"scalefactor" : 0.02,
	"crossaxis_sensitivity" : 0.02,
	"b" : 40e3,  # initialize with np.random and then don't change ever
	"cov" : 0.0005
}

sun_params = {
	"scalefactor" : 0,
	"crossaxis_sensitivity" : 0.02,
	"b" : 1.0, # initialize with np.random and then don't change ever
	"cov" : 0.0005
}

mag_order = 10
