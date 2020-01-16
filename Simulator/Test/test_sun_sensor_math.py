# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 17:22:15 2019

"""
import pytest
import os, sys, inspect
import sun_sensor_math as sensors
import astropy_sun_position as sun
import numpy as np

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)

def test_sense2vector():

    meas = [1, 2, 0, 0, 0, 0]
    r_sat = [8000, 0, 0]
    q = [1, 0, 0, 0] # identity quaternion

    vec = sensors.sense2vector(meas, r_sat, q, albedo = False)

    assert np.testing.assert_allclose(vec, [-1, 0, 0], atol = 10e-5)

    ######### use full information to compute measurements from known sun vector and then go back and see if they match. ####
    # r_sat = [6371., 2123., 1061.]
    # q = quat(np.array([1,2,3,4]))

    # mjd = 54000
    # r_Earth2Sun = sun.sun_position_ECI(mjd)

    # vec = sensors.sense2vector(meas, r_sat, q)

    # assert




def test_deltas2measure():
    deltas = [1,-1,1]
    b = all(LA.deltas2measure(deltas) == [1,0,0,1,1,0])
    assert b == True

    deltas = [-1,-1,1]
    b = all(LA.deltas2measure(deltas) == [0,1,0,1,1,0])
    assert b == True

    deltas = [-100,1,1]
    b = all(LA.deltas2measure(deltas) == [0,100,1,0,1,0])
    assert b == True

    deltas = [-100,1]
    b = all(LA.deltas2measure(deltas) == [0,100,1,0])
    assert b == True

    deltas = [-100]
    b = all(LA.deltas2measure(deltas) == [0,100])
    assert b == True
