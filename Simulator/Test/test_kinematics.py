import pytest
import os, sys, inspect, pdb
import numpy as np

# add current folder to the path
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)

from kinematics import calc_q_dot
from conversions import L, quat
tol = 1e-6


def test_calc_q_dot():
	q = quat(np.random.rand(4))
	w = np.random.rand(3)
	test1 = .5 * L(q) @ np.concatenate((np.array([0]), w))
	test2 = calc_q_dot(q, w)

	np.testing.assert_allclose(test1, test2, atol=tol)