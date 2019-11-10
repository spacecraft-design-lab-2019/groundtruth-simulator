# -*- coding: utf-8 -*-

"""
Verify against SGP4
"""

import numpy as np
from sim_config import *
from simulation_step import simulation_step
from propagate_step import sgp4_step


# Initialize / Setup Workspace
tspan = np.array([0, 8640])    # [sec]
T = np.arange(0, tspan[1]+tstep, tstep)


# Initial State Vector
r_i, v_i = sgp4_step(line1, line2, tstart)
state_i = np.r_[r_i, q_i, v_i, w_i]
state_history = np.zeros((np.shape(T)[0], np.shape(state_i)[0]))
state_history_sgp4 = np.zeros((np.shape(T)[0], 6))
state_history[0, :] = state_i
state_history_sgp4[0, :] = np.r_[r_i, v_i]
sim_state = {'state': state_i, 't': tstart}



# Propagate
for i, elapsed_t in enumerate(T[0:-1]):

	# Simulator
	sensors, sim_state = simulation_step(np.zeros((3,1)), sim_state)
	state_history[i+1, :] = sim_state['state']

	# SGP4
	state_history_sgp4[i+1, :] = np.array(sgp4_step(line1, line2, sim_state['t']))


# Plot
plt.figure()
plt.plot(T, abs(state_history[:, 0:3] - state_history_sgp4[:, 0:3]))
plt.plot(T, abs(state_history[:, 7:10] - state_history_sgp4[:, 3:6]), hold='on')
plt.xlabel('time')
plt.ylabel('error')
plt.suptitle('Error against SGP4')
