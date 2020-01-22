import pytest
import os, sys, inspect, pdb
import numpy as np

# add current folder to the path
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)

import simulator
from constants import SpacecraftStructure, Environment
import julian
import sim_config as config


def test_euler():
	# test euler equations + magnetic dipole for 60 second run against Andrew's sim

	# STARTING CONDITIONS
	r_i = np.array([0.5001945144886031, 1.0083208747495578, 0.0]) * 6378.1378366 # km
	v_i = np.array([-3.588574367865696, 1.780172620235258, 6.287934046067301]) # km/s
	q_i = np.array([1, 0, 0, 0])
	w_i = np.array([.1, .1, .1]) # rad/s
	mjd_start = 58865.5
	I = np.array([[.15, .015, .015],[.015, .15, .015],[.015, .015, .15]])
	mass = 10.0 # kg
	L_cmd = np.ones(3) # [Am^2 or Nm/T]

	# SETUP SIM - note we have to re-do work from the __init__ function to avoid editing sim_config just for a unit_test
	sim = simulator.Simulator(config)
	sim.state = np.r_[r_i, q_i, v_i, w_i]
	sim.t = julian.from_jd(mjd_start, fmt='mjd')
	sim.MJD = mjd_start
	sim.tstep = .1
	sim.environment = Environment(sim.t)
	sim.structure = SpacecraftStructure(I, mass=mass)

	# PREALLOCATE MEMORY
	tspan = np.array([0, 60])    # [sec]
	T = np.arange(0, tspan[1]+sim.tstep, sim.tstep)
	state_history = np.zeros((np.shape(T)[0], np.shape(sim.state)[0]))
	state_history[0, :] = sim.state

	# RUN SIM
	for i, elapsed_t in enumerate(T[0:-1]):
		sim.step(config.tstep, L_cmd)
		state_history[i+1, :] = sim.state

	# ANSWER - from Andrew's sim
	w_ans = np.array([0.10192649869210156, 0.09085024017514931, 0.10722326113274933]) # rad/s
	q_ans = np.array([0.47006180334049474, -0.5230023039391511, -0.4832166552830285, -0.5215478455215423])
	r_ans = np.array([0.46544750473390517, 1.0231138587560542, 0.05917931343363079]) * 6378.1378366 # km
	v_ans = np.array([-3.7881391987173676, 1.3603452712714288, 6.275703727336128]) # km/s

	w_sim = state_history[-1, 10:13]
	q_sim = state_history[-1, 3:7]
	r_sim = state_history[-1, 0:3]
	v_sim = state_history[-1, 7:10]

	print(w_ans)
	print(w_sim)

	print(q_sim)
	print(q_ans)

	print(r_sim)
	print(r_ans)

	print(v_sim)
	print(v_ans)


	np.testing.assert_allclose(w_sim, w_ans, atol=1e-2)
	np.testing.assert_allclose(q_sim, q_ans, atol=1e-1)
	np.testing.assert_allclose(r_sim, r_ans, atol=1e0)
	np.testing.assert_allclose(v_sim, v_ans, atol=1e-2)