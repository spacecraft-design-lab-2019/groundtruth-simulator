# -*- coding: utf-8 -*-

import numpy as np
import sim_config as config
from constants import *
from propagate_step import *


def simulation_step(cmd, sim_prev=None):
	"""
	Function: simulation_step
		Propagates dynamics & models sensors for single step

	Inputs:
		cmd:		commanded magnetic dipole from controller
		sim_prev: 	previous state of simulation

	Outputs:
		sensors: 	spoofed sensor measurements
		sim_new:	new state of simulation
	"""

	#------------------ Initialize/Setup Workspace ------------------	
	# if we are at the first iteration
	if sim_prev == None:
		r_i, v_i = sgp4_step(line1, line2, tstart)
		state_i = np.r_[r_i, config.q_i, v_i, config.w_i]
		t_i = config.tstart

	# get previous state (groundtruth)
	else:
		state_i = sim_prev.state
		t_i = sim_prev.t

	# environment + spacecraft classes
	world_i = Environment(t_i)
	spacecraft = SpacecraftStructure(config.I)


	#------------------------ Propagate Dynamics --------------------
	update_f = lambda t, state: calc_statedot(t, state, cmd, world_i, spacecraft)
	state = rk4_step(update_f, t_i, state_i, config.tstep)


	#------------------------ Spoof Sensors -------------------------
	# TO-DO: use sensor classes to spoof sensors based on updated world and state
	world = Environment(t_i + datetime.timedelta(seconds=config.tstep))
	sensors = []

	#------------------------ Export Data -------------------------
	# TO-DO: output desired variables to text file for later plotting/analysis
	sim_new = []

	return sensors, sim_new