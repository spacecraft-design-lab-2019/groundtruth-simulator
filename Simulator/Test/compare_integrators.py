# -*- coding: utf-8 -*-

import numpy as np
import time
import sys, os
import sim_config as config
from simulator import Simulator


#-----------------Configuration / Parameters--------------------
tspan = np.array([0, 100])    # [sec]


#----------------Initialize / Setup Workspace------------------
# setup sim
sim = Simulator(config)

# preallocate memory
T = np.arange(0, tspan[1]+config.tstep, config.tstep)

#---------------------Propagate---------------------------
t = time.time()
for i, elapsed_t in enumerate(T[0:-1]):
	sim.step()
elapsed = time.time() - t

print(elapsed)