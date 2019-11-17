# -*- coding: utf-8 -*-

import numpy as np
import conversions as conv
from propagate_step import *
from constants import SpacecraftStructure, Environment
from sensormodels import SpacecraftSensors


class Simulator():
	"""
	A class to initialize and run the simulator.
	"""
	def __init__(self, config):

		# initialize classes
		self.structure = SpacecraftStructure(config.I, mass=config.mass)
		self.environment = Environment(config.tstart)
		self.sensors = SpacecraftSensors(config.mag_params, config.gyro_params, config.sun_params);

		# initial state
		r_i, v_i = sgp4_step(config.line1, config.line2, config.tstart)
		self.state = np.r_[r_i, config.q_i, v_i, config.w_i]
		self.t = config.tstart


	def step(self, cmd, tstep):
		"""
		Function: step
			Propagates dynamics & models sensors for single step

		Inputs:
			cmd:	commanded magnetic dipole from controller
		Outputs:
			meas: 	spoofed sensor measurements

		Comments:
			-Simulation state format: [ r [km], q [], v [km/s], w [rad/s] ] (13x1 np array)
			-r & v are in ECI coordinates, q is the transformation from body to ECI, and w is given in the body frame

		"""
		#------------------------ Propagate Dynamics --------------------
		update_f = lambda t, state: calc_statedot(t, state, cmd, self.structure, self.environment)
		self.state = rk4_step(update_f, self.t, self.state, tstep)
		self.t = self.t + datetime.timedelta(seconds=tstep)


		#------------------------ Calculate Environment -------------------
		self.environment.update(self.t)

		B_ECI = self.environment.magfield_lookup(self.state[0:3])
		B_body = conv.quatrot(self.state[3:7], B_ECI)

		# TO-DO: update gyro bias
		#self.environment.sensors.gyroscope.errormodel.b += ?
		w_body = conv.quatrot(self.state[3:7], self.state[10:13])


		#------------------------ Spoof Sensors -------------------------
		B_body_noise = self.sensors.magnetometer.measure(B_body)
		w_body_noise = self.sensors.gyroscope.measure(w_body)

		meas = np.r_[B_body_noise, w_body_noise]

		#------------------------ Export Data -------------------------
		# TO-DO: output desired variables to text file for later plotting/analysis

		return meas