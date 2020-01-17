# -*- coding: utf-8 -*-

import numpy as np
from datetime import datetime
import sys
from math import sqrt
sys.path.append('/home/eleboeuf/Documents/GNC')
sys.path.append('./GNC')
import time_functions_cpp
#------------------- Configuration Parameters -------------------

# Seed Initial Position/Velocity with TLE - BEESAT-1
# (optional) - can instead replace this with r_i, v_i as np.array(3)
line1 = ('1 35933U 09051C   19315.45643387  .00000096  00000-0  32767-4 0  9991')
line2 = ('2 35933  98.6009 127.6424 0006914  92.0098 268.1890 14.56411486538102')

# Simulation Parameters
tstart = datetime(2019, 12, 30, 00, 00, 00)
tstep = .1                     # [sec] - 1 Hz
MJDstart = time_functions_cpp.date2MJD(12, 30, 2019, 0, 0, 0)
# Initial Spacecraft Attitude
q_i = np.array([sqrt(4.0)/4.0, sqrt(4.0)/4.0, sqrt(4.0)/4.0, sqrt(4.0)/4.0])    # quaternion
w_i = np.array([.03, .03, .03])   # radians/sec

# Spacecraft Properties
I = np.array([[.3,0,0],[0,.4,0],[0,0,.5]])
mass = 1.0 # kg

# Sensor Parameters
gyro_params = {
	"scaleF" : 0.002,
	"caSense" : 0.02,
	"b" : np.sqrt(.0022), # will be mulitplied to np.random.rand(3)
	"cov" : 0.000000000694444
}

mag_params = {
	"scaleF" : 0.02,
	"caSense" : 0.02,
	"b" : 40e3,  # initialize with np.random and then don't change ever
	"cov" : 0.0005
}

sun_params = {
	"scaleF" : 0,
	"caSense" : 0.02,
	"b" : 1.0, # initialize with np.random and then don't change ever
	"cov" : 0.0005
}

mag_order = 10