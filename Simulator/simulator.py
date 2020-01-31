# -*- coding: utf-8 -*-

import numpy as np
import datetime, pdb
import conversions as conv
from scipy.integrate import solve_ivp
from propagate_step import sgp4_step, rk4_step, calc_statedot
from constants import SpacecraftStructure, Environment
from sensormodels import SpacecraftSensors

class Simulator():
	"""
	A class to initialize and run the simulator.
	"""
	def __init__(self, config):

		# initialize classes
		self.structure = SpacecraftStructure(config.I, mass=config.mass)
		self.environment = Environment(config.mjd_start)
		self.sensors = SpacecraftSensors(config.mag_params, config.gyro_params, config.sun_params);

		# simulation time [sec]
		self.t = 0
		self.tstep = config.tstep

		# initial state
		r_i, v_i = sgp4_step(config.line1, config.line2, config.mjd_start)
		self.state = np.r_[r_i, config.q_i, v_i, config.w_i]

		# magnetic field order (for IGRF)
		self.mag_order = config.mag_order


	def step(self, cmd=np.zeros(3)):
		"""
		Function: step
			Propagates dynamics & models sensors for single step

		Inputs:
			cmd:	commanded magnetic moment, [Am^2 or Nm/T]
		Outputs:
			meas: 	spoofed sensor measurements

		Comments:
			-Simulation state format: [ r [km], q [], v [km/s], w [rad/s] ] (13x1 np array)
			-Command format: Torque [Nm]
			-r & v are in ECI coordinates, q is the transformation from body to ECI, and w is given in the body frame

		"""
		#------------------------ Propagate Dynamics --------------------
		update_f = lambda t, state: calc_statedot(t, state, cmd, self.structure, self.environment, self.mag_order)
		# sol = solve_ivp(update_f, (self.t, self.t+self.tstep), self.state)
		# self.t = sol.t[-1]
		# self.state = sol.y[:,-1]
		self.t, self.state = rk4_step(update_f, self.t, self.state, self.tstep)
		
		self.state[3:7] = self.state[3:7] / np.linalg.norm(self.state[3:7]) # normalize the quaternion vector

		#------------------------ Calculate Environment -------------------
		B_ECI = self.environment.magfield_lookup(self.state[0:3], self.mag_order) # Earth's magnetic field isn't fixed in ECI space, it's fixed in ECEF space!!!!
		B_body = conv.quatrot(conv.conj(self.state[3:7]), B_ECI)

		S_ECI = self.environment.sunVector(self.state[0:3])		
		S_body = conv.quatrot(conv.conj(self.state[3:7]), S_ECI)


		#------------------------ Spoof Sensors -------------------------
		# Actuate based on truth for now until magnetometer bias estimation, TRIAD, and MEKF have been implemented and tested
		# B_body_noise = B_body
		S_body_noise = S_body
		# S_body_noise = self.sensors.sunsensor.measure(S_body)
		w_body_noise = self.state[10:13]
		B_body_noise = self.sensors.magnetometer.measure(B_body)
		# w_body_noise = self.sensors.gyroscope.measure(self.state[10:13])

		meas = np.r_[B_body_noise, w_body_noise, S_body_noise]

		#------------------------ Export Data -------------------------
		# TO-DO: output desired variables to text file for later plotting/analysis
		self.debug_output = [B_ECI, B_body]

		return meas