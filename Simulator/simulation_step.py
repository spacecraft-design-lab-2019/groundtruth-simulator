# -*- coding: utf-8 -*-

import numpy as np
import sim_config as config
import conversions as conv
from constants import *
from propagate_step import *
import sensors as sense


def simulation_step(cmd, sim_prev=None):
	"""
	Function: simulation_step
		Propagates dynamics & models sensors for single step

	Inputs:
		cmd:		commanded magnetic dipole from controller
		sim_prev: 	previous state of simulation

	Outputs:
		sensors: 	spoofed sensor measurements
		sim_new:	new state of simulation, dictionary.

	Comments:
		-Simulation state format: [ r [km], q [], v [km/s], w [rad/s] ] (13x1 np array)
		-r & v are in ECI coordinates, q is the transformation from body to ECI, and w is given in the body frame
		-
	"""

	#------------------ Initialize/Setup Workspace ------------------	
	# if we are at the first iteration
	if sim_prev == None:
		r_i, v_i = sgp4_step(config.line1, config.line2, config.tstart)
		state_i = np.r_[r_i, config.q_i, v_i, config.w_i]
		t_i = config.tstart

	# get previous state (groundtruth)
	else:
		state_i = sim_prev['state']
		t_i = sim_prev['t']

	# spacecraft class
	spacecraft = SpacecraftStructure(config.I, mass=config.mass)


	#------------------------ Propagate Dynamics --------------------
	update_f = lambda t, state: calc_statedot(t, state, cmd, spacecraft)
	state = rk4_step(update_f, t_i, state_i, config.tstep)
	t = t_i + datetime.timedelta(seconds=config.tstep)

	#------------------------ Extract true values -------------------
	# TODO: change so that we don't create a new environment object every step. Add "increment_time" method and pass environment
	# into/out of this simulation step function
	world = Environment(t)

	B_NED = world.magfield_lookup(state[0:3])

	# TODO: Add conversion from NED to ECI


	B_body = conv.quatrot(state[3:7], B_ECI)
	w_body = conv.quatrot(state[3:7], state[10:13])

	#------------------------ Spoof Sensors -------------------------
	# TODO: use sensor classes to spoof sensors based on updated world and state

	# Uncomment 2nd line to include sensor models
	B_body_noise = B_body
	# B_body_noise = sense.magnetometerModel(B_body)

	# Uncomment 2nd line to include sensor models
	w_body_noise = w_body
	# w_body_noise = sense.gyroModel(w_body)
	sensors = np.r_[B_body_noise, w_body_noise]

	#------------------------ Export Data -------------------------
	# TO-DO: output desired variables to text file for later plotting/analysis
	sim_new = {'state': state, 't': t}

	return sensors, sim_new