import pytest
import os, sys, inspect
import numpy as np

# add current folder to the path
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)

import conversions as conv
tol = 1e-6

def test_ECI_to_ECEF():
    # Case 1
    GMST = 0
    x = np.random.randn(3)
    np.testing.assert_allclose(x, conv.ECI_to_ECEF(x, GMST), atol=tol)

    # Case 2
    GMST = np.pi/2
    x = np.array([0, 1, 0])
    y = np.array([1, 0, 0])
    np.testing.assert_allclose(y, conv.ECI_to_ECEF(x, GMST), atol=tol)

def test_ECEF_to_ECI():
    # Case 1
    GMST = 0
    x = np.random.randn(3)
    np.testing.assert_allclose(x, conv.ECEF_to_ECI(x, GMST), atol=tol)

    # Case 2
    GMST = np.pi/2
    x = np.array([0, 1, 0])
    y = np.array([1, 0, 0])
    np.testing.assert_allclose(x, conv.ECEF_to_ECI(y, GMST), atol=tol)

def test_ECEF_to_LLA():
    x = np.array([1, 0, 0])
    lla = [0, 0, -6378.1378366]
    np.testing.assert_allclose(lla, conv.ECEF_to_LLA(x, 6378.1378366), atol=tol)