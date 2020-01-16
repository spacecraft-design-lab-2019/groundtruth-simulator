# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 17:22:15 2019

"""
import pytest
import os, sys, inspect
import sun_sensor_math as sensors
import astropy_sun_position as sun
import numpy as np
import conversions as conv

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)

def test_sense2vector():

    meas = [1., 2., 0., 0., 0., 0.]
    r_sat = [8000, 0, 0]
    q = [1, 0, 0, 0] # identity quaternion

    vec = sensors.sense2vector(meas, r_sat, q, albedo = False)

    np.testing.assert_allclose(vec, [-1., 0., 0.], atol = 10e-5)

    ######### use full information to compute measurements from known sun vector and then go back and see if they match. ####
    # r_sat = [6371., 2123., 1061.]
    # q = quat(np.array([1,2,3,4]))

    mjd = 54000
    r_Earth2Sun = sun.sun_position_ECI(mjd)
    r_sat = [8000,0,0]
    q = [1,1,2,1] #eci to body quaternion
    sat2sun_input = conv.quatrot(q,sensors.add(r_Earth2Sun,r_sat)) #in body frame
    meas = sensors.vector2sense(sat2sun_input, r_sat, q, albedo = False)
    sat2sun_out = sensors.sense2vector(meas, r_sat, q, albedo = False) #in body frame
    out_r_Earth2Sun = sensors.sub(conv.quatrot(conv.conj(q),sat2sun_out),r_sat) #output in eci
    
    np.testing.assert_allclose(sat2sun_out, sensors.normalize(sat2sun_input), atol = 10e-5)
#    np.testing.assert_allclose(out_r_Earth2Sun, r_Earth2Sun, atol = 10e-5)

    in_meas = [1., 0., 0., 0., 0., 0.]
    r_sat = [8000, 0, 0]
    q = [1, 0, 0, 0] # identity quaternion
    sat2sun = sensors.sense2vector(in_meas, r_sat, q, albedo = False)
    out_meas = sensors.vector2sense(sat2sun, r_sat, conv.conj(q), albedo = False)
    
    np.testing.assert_allclose(out_meas, in_meas, atol = 10e-5)
    
    
    in_meas = [2, 1, 100, 100, 0, 0]
    sat2sun = sensors.sense2vector(in_meas, r_sat, q, albedo = False)
    out_meas = sensors.vector2sense(sat2sun, r_sat, conv.conj(q), albedo = False)
    
    np.testing.assert_allclose(out_meas, [1,0,0,0,0,0], atol = 10e-5)


def test_deltas2measure():
    deltas = [1,-1,1]
    b = all(sensors.deltas2measure(deltas) == [1,0,0,1,1,0])
    assert b == True

    deltas = [-1,-1,1]
    b = all(sensors.deltas2measure(deltas) == [0,1,0,1,1,0])
    assert b == True

    deltas = [-100,1,1]
    b = all(sensors.deltas2measure(deltas) == [0,100,1,0,1,0])
    assert b == True

    deltas = [-100,1]
    b = all(sensors.deltas2measure(deltas) == [0,100,1,0])
    assert b == True

    deltas = [-100]
    b = all(sensors.deltas2measure(deltas) == [0,100])
    assert b == True

