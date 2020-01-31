import pytest
import os, sys, inspect, pdb
import numpy as np

# add current folder to the path
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)
sys.path.append('./GNC')  # this line is specifically for the GitHub validation check that runs on their servers

import simulator
from constants import SpacecraftStructure, Environment
import julian
import sim_config as config


def test_dynamics():
	# test orbit + attitude dynamics for 60 second run against Andrew's sim

	# STARTING CONDITIONS
	r_i = np.array([0.5006827813414019, 1.009305151157764, 0.0]) * 6378.1378366 # km
	v_i = np.array([-3.5868241485959342, 1.7793043945545146, 6.284867295261474]) # km/s
	q_i = np.array([1, 0, 0, 0])
	w_i = np.array([.1, .1, .1]) # rad/s
	mjd_start = 58865.5
	I = np.array([[.15, .015, .015],[.015, .15, .015],[.015, .015, .15]])
	mass = 10.0 # kg
	L_cmd = np.ones(3) # [Am^2 or Nm/T]

	# SETUP SIM - note we have to re-do work from the __init__ function to avoid editing sim_config just for a unit_test
	sim = simulator.Simulator(config)
	sim.state = np.r_[r_i, q_i, v_i, w_i, config.T_i]
	sim.t = julian.from_jd(mjd_start, fmt='mjd')
	sim.MJD = mjd_start
	sim.tstep = .1
	sim.environment = Environment(sim.t)
	sim.structure = SpacecraftStructure(I, mass, config.thermal_properties)

	# PREALLOCATE MEMORY
	tspan = np.array([0, 60])    # [sec]
	T = np.arange(0, tspan[1]+sim.tstep, sim.tstep)
	state_history = np.zeros((np.shape(T)[0], np.shape(sim.state)[0]))
	state_history[0, :] = sim.state

	bECI_history = np.zeros((np.shape(T)[0], 3))

	# RUN SIM
	for i, elapsed_t in enumerate(T[0:-1]):
		sim.step(config.tstep, L_cmd)
		state_history[i+1, :] = sim.state
		bECI_history[i+1, :] = sim.debug_output[0]/1e9

	# ANSWER - from Andrew's sim
	w_ans = np.array([0.103473377474831, 0.10351893953023031, 0.0930076829949389]) # rad/s
	q_ans = np.array([0.46758555857559125, -0.5129972392025514, -0.5249064223617074, -0.49261630681699486])
	r_ans = np.array([0.4659942472475783, 1.024081911258644, 0.05908576934024495]) * 6378.1378366 # km
	v_ans = np.array([-3.785581159645769, 1.3612227557429792, 6.272718772389778]) # km/s

	w_sim = state_history[-1, 10:13]
	q_sim = state_history[-1, 3:7]
	r_sim = state_history[-1, 0:3]
	v_sim = state_history[-1, 7:10]

	np.testing.assert_allclose(w_sim, w_ans, atol=1e-3)
	np.testing.assert_allclose(q_sim, q_ans, atol=1e-3)
	np.testing.assert_allclose(r_sim, r_ans, atol=1e-1)
	np.testing.assert_allclose(v_sim, v_ans, atol=1e-2)

	# note that andrew's sim did not account for aero drag/torque or gravity gradient torques
	# so, these tolerances are reasonably chosen