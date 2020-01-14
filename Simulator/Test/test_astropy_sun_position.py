
import pytest
import os, sys, inspect

# Add current folder to the path
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)

import astropy_sun_position as asp
tol = 1e-6


def test_approx_sun_position_ECI():
	assert True

def test_sun_position_ECI():
	assert True

def test_GCRS_to_ECI():
	assert True

