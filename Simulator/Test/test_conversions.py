import pytest
import os, sys, inspect
import numpy as np

# add current folder to the path
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)

import conversions as conv
tol = 1e-6

#--------------------Coordinate Frames--------------------------

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
    # Case 1
    x = np.array([1, 0, 0])
    lla = [0, 0, 0]
    np.testing.assert_allclose(lla, conv.ECEF_to_LLA(x, 1), atol=tol)

    # Case 2
    x = np.array([0, 0, 1])
    lla = [np.pi/2, 0, 0]
    np.testing.assert_allclose(lla, conv.ECEF_to_LLA(x, 1), atol=tol)

    # Case 3
    x = np.array([0, 1, 0])
    lla = [0, np.pi/2, 0]
    np.testing.assert_allclose(lla, conv.ECEF_to_LLA(x, 1), atol=tol)

def test_NED_to_ECI():
    vec_NED = np.array([0, 0, 1])
    vec_ECI = np.array([-1, 0, 0])
    np.testing.assert_allclose(vec_ECI, conv.NED_to_ECI(vec_NED, 0, 0, 0), atol=tol)


#--------------------Quaternions-----------------------------

def test_quatrot():
    # Verified against Matlab quatrotate()
    vec = np.array([-0.7734, -1.7788, -0.9937])
    quat = conv.conj(np.array([1, 0, 1, 0]))
    ans = np.array([0.9937, -1.7788, -0.7734])
    np.testing.assert_allclose(ans, conv.quatrot(quat/np.linalg.norm(quat), vec), atol=tol)

    vec = np.array([1, 1, 0])
    quat = np.array([1, 0, 1, 0])
    ans = np.array([0, 1, -1])
    np.testing.assert_allclose(ans, conv.quatrot(quat/np.linalg.norm(quat), vec), atol=tol)

    vec = np.array([1, 1, 0])
    quat = np.array([1, 0, 1, 0])
    ans = np.array([0, 1, -1])
    np.testing.assert_allclose(ans, conv.quatrot(quat/np.linalg.norm(quat), vec), atol=tol)
    

#--------------------Miscellaneous-----------------------------

def test_skew():
    x = np.array([1, 0, 0])
    y = np.array([0, 1, 0])
    z = np.array([0, 0, 1])
    zero = np.zeros(3)
    np.testing.assert_allclose(z, conv.skew(x)@y, atol=tol)
    np.testing.assert_allclose(zero, conv.skew(x)@x, atol=tol)

def test_mjd_2_GMST():
    assert True

def test_unit():
    np.testing.assert_allclose(np.array([1,0,0]), conv.unit('x'), atol=tol)
    np.testing.assert_allclose(np.array([0,1,0]), conv.unit('y'), atol=tol)
    np.testing.assert_allclose(np.array([0,0,1]), conv.unit('z'), atol=tol)