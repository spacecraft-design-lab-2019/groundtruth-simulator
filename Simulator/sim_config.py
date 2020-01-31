# -*- coding: utf-8 -*-
import os, sys
import julian
import numpy as np
from datetime import datetime
from math import sqrt

#------------------- Configuration Parameters -------------------

# Simulation Parameters
tstart = julian.to_jd(datetime(2020, 1, 30, 00, 00, 00), fmt='mjd')
tstep = 0.1 / (60*60*24)  # 10 Hz


# Seed Initial Position/Velocity with TLE - BEESAT-1
# (optional) - can instead replace this with r_i, v_i as np.array(3)
line1 = ('1 35933U 09051C   19315.45643387  .00000096  00000-0  32767-4 0  9991')
line2 = ('2 35933  98.6009 127.6424 0006914  92.0098 268.1890 14.56411486538102')


# Initial Spacecraft Attitude
q_i = np.array([sqrt(4.0)/4.0, sqrt(4.0)/4.0, sqrt(4.0)/4.0, sqrt(4.0)/4.0])    # quaternion
w_i = np.array([.03, .03, .03])   # radians/sec


# Spacecraft Properties
I = np.array([[.3,0,0],[0,.4,0],[0,0,.5]])
mass = 1.0 # kg


# Sensor Parameters
mag_order = 10
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