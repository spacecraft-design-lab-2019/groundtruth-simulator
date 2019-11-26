# -*- coding: utf-8 -*-

import numpy as np
import conversions as conv
from propagate_step import *
from constants import SpacecraftStructure, Environment
from sensormodels import SpacecraftSensors
import sys
sys.path.append('/home/eleboeuf/Documents/GNC')
import sun_utils_cpp

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
		self.MJD = config.MJDstart
		self.tstep = config.tstep


	def step(self, tstep, cmd=np.zeros(3)):
		"""
		Function: step
			Propagates dynamics & models sensors for single step

		Inputs:
			cmd:	commanded torque, [Nm]
		Outputs:
			meas: 	spoofed sensor measurements

		Comments:
			-Simulation state format: [ r [km], q [], v [km/s], w [rad/s] ] (13x1 np array)
			-Command format: Torque [Nm]
			-r & v are in ECI coordinates, q is the transformation from body to ECI, and w is given in the body frame

		"""
		#------------------------ Propagate Dynamics --------------------
		update_f = lambda t, state: calc_statedot(t, state, cmd, self.structure, self.environment)
		self.state = rk4_step(update_f, self.t, self.state, tstep)
		self.t = self.t + datetime.timedelta(seconds=tstep)
		self.MJD = self.MJD + self.tstep / 24 / 3600

		#------------------------ Calculate Environment -------------------
		self.environment.update(self.t)
		self.sensors.gyroscope.update_bias()

		B_ECI = self.environment.magfield_lookup(self.state[0:3])
		B_body = conv.quatrot(self.state[3:7], B_ECI)

		# S_ECI = sun_utils_cpp.sun_position(self.MJD) 
		S_ECI = sun_utils_cpp.sat_sun_vect(self.state[0:3], self.MJD) 
		S_ECI = S_ECI / np.linalg.norm(S_ECI)
		S_body = conv.quatrot(self.state[3:7], S_ECI)


		#------------------------ Spoof Sensors -------------------------
		# Actuate based on truth for now until magnetometer bias estimation, TRIAD, and MEKF have been implemented and tested
		B_body_noise = B_body
		S_body_noise = S_body
		w_body_noise = self.state[10:13]
		# B_body_noise = self.sensors.magnetometer.measure(B_body)
		# w_body_noise = self.sensors.gyroscope.measure(self.state[10:13])

		meas = np.r_[B_body_noise, w_body_noise, S_body_noise]

		#------------------------ Export Data -------------------------
		# TO-DO: output desired variables to text file for later plotting/analysis
		self.debug_output = [B_ECI, B_body]

		return meas