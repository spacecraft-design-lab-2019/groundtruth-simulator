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
    
    meas = [1., 2., 0., 0., 0., 0.]
    r_sat = [8000, 0, 0]
    q = [1, 0, 0, 0] # identity quaternion

    vec = sensors.sense2vector(meas, r_sat, q, albedo = True)

    np.testing.assert_allclose(vec, [-1., 0., 0.], atol = 10e-5)
    
    meas = [59., 60., 2., 1., 1., 0.]
    r_sat = [8000, 0, 0]
    q = [1, 0, 0, 0] # identity quaternion

    vec = sensors.sense2vector(meas, r_sat, q, albedo = False)

    np.testing.assert_allclose(vec, sensors.normalize([-1., 1., 1.]), atol = 10e-5)

    mjd = 54000
    r_Earth2Sun = sun.sun_position_ECI(mjd)
    r_sat = [8000,0,0]
    q = [1,1,2,1] #eci to body quaternion
    sat2sun_input = conv.quatrot(q,sensors.add(r_Earth2Sun,r_sat)) #in body frame
    meas = sensors.vector2sense(sat2sun_input, r_sat, q, albedo = False)
    sat2sun_out = sensors.sense2vector(meas, r_sat, q, albedo = False) #in body frame
    
    np.testing.assert_allclose(sat2sun_out, sensors.normalize(sat2sun_input), atol = 10e-5)
    
    mjd = 54104
    r_Earth2Sun = sun.sun_position_ECI(mjd)
    r_sat = [7543,201,345]
    q = [1,1,0,1] #eci to body quaternion
    sat2sun_input = conv.quatrot(q,sensors.add(r_Earth2Sun,r_sat)) #in body frame
    meas = sensors.vector2sense(sat2sun_input, r_sat, q, albedo = False)
    sat2sun_out = sensors.sense2vector(meas, r_sat, q, albedo = False) #in body frame
    
    np.testing.assert_allclose(sat2sun_out, sensors.normalize(sat2sun_input), atol = 10e-5)
    
    mjd = 52124
    r_Earth2Sun = sun.sun_position_ECI(mjd)
    r_sat = [543,201,7345]
    q = [1,.1,.3,1] #eci to body quaternion
    sat2sun_input = conv.quatrot(q,sensors.add(r_Earth2Sun,r_sat)) #in body frame
    meas = sensors.vector2sense(sat2sun_input, r_sat, q, albedo = False)
    sat2sun_out = sensors.sense2vector(meas, r_sat, q, albedo = False) #in body frame
    
    np.testing.assert_allclose(sat2sun_out, sensors.normalize(sat2sun_input), atol = 10e-5)
    
    mjd = 54134.1
    r_Earth2Sun = sun.sun_position_ECI(mjd)
    r_sat = [1073,7000,345]
    q = conv.quat([.2,.8,.4,.1]) #eci to body quaternion
    sat2sun_input = conv.quatrot(q,sensors.add(r_Earth2Sun,r_sat)) #in body frame
    meas = sensors.vector2sense(sat2sun_input, r_sat, q, albedo = True)
    sat2sun_out = sensors.sense2vector(meas, r_sat, q, albedo = True) #in body frame
    
    np.testing.assert_allclose(sat2sun_out, sensors.normalize(sat2sun_input), atol = 15e-3) #acceptable deviation
    
    r_Earth2Sun = [1.5018e08,0,0]
    r_sat = [0,7000,0]
    q = conv.quat([1,.2,.8,.2]) #eci to body quaternion
    sat2sun_input = conv.quatrot(q,sensors.add(r_Earth2Sun,r_sat)) #in body frame
    meas = sensors.vector2sense(sat2sun_input, r_sat, q, albedo = True)
    sat2sun_out = sensors.sense2vector(meas, r_sat, q, albedo = True) #in body frame
    
    np.testing.assert_allclose(sat2sun_out, sensors.normalize(sat2sun_input), atol = 50e-4)
    
    r_Earth2Sun = [0,1.5018e08,0]
    r_sat = [0,7000,0]
    q = conv.quat([1,.2,.8,.2]) #eci to body quaternion
    sat2sun_input = conv.quatrot(q,sensors.add(r_Earth2Sun,r_sat)) #in body frame
    meas = sensors.vector2sense(sat2sun_input, r_sat, q, albedo = True)
    sat2sun_out = sensors.sense2vector(meas, r_sat, q, albedo = True) #in body frame
    
    np.testing.assert_allclose(sat2sun_out, sensors.normalize(sat2sun_input), atol = 1e-5)
    
    mjd = 54134.9
    r_Earth2Sun = sun.sun_position_ECI(mjd)
    r_sat = [523,6201,345]
    q = conv.quat([-.6,1,3,1]) #eci to body quaternion
    sat2sun_input = conv.quatrot(q,sensors.add(r_Earth2Sun,r_sat)) #in body frame
    meas = sensors.vector2sense(sat2sun_input, r_sat, q, albedo = True)
    sat2sun_out = sensors.sense2vector(meas, r_sat, q, albedo = True) #in body frame
    
    np.testing.assert_allclose(sat2sun_out, sensors.normalize(sat2sun_input), atol = 15e-3)
    
    r_Earth2Sun = [1.5018e08,0,0]
    r_sat = [0,7000,0]
    q = conv.quat([1,1,1,1]) #eci to body quaternion
    sat2sun_input = conv.quatrot(q,sensors.add(r_Earth2Sun,r_sat)) #in body frame
    meas = sensors.vector2sense(sat2sun_input, r_sat, q, albedo = True)
    sat2sun_out = sensors.sense2vector(meas, r_sat, q, albedo = True) #in body frame
    
    np.testing.assert_allclose(sat2sun_out, sensors.normalize(sat2sun_input), atol = 5e-3)
    
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

def test_isEclipse():
    Re = 6371
    
    r_sat = [8000,0,0]
    r_Earth2Sun = [-1.50147817e8,  2.97758769e6, -2.55999125e4]
    assert sensors.isEclipse(r_sat, r_Earth2Sun, Re) == True
    
    r_sat = [8000,0,0]
    r_Earth2Sun = [-2e6,  1.50147817e8,  2.97758769e6]
    assert sensors.isEclipse(r_sat, r_Earth2Sun, Re) == False
    
    r_sat = [8000,0,0]
    r_Earth2Sun = [-95e6,  1.50147817e8,  2.97758769e6]
    assert sensors.isEclipse(r_sat, r_Earth2Sun, Re) == False
    
    r_sat = [8000,0,0]
    r_Earth2Sun = [-120e6,  1.50147817e8,  2e6]
    assert sensors.isEclipse(r_sat, r_Earth2Sun, Re) == True
    
    r_sat = [-8000,0,0]
    r_Earth2Sun = [120e6,  1.50147817e8,  2e6]
    assert sensors.isEclipse(r_sat, r_Earth2Sun, Re) == True
    
    r_sat = [600,-7800,-500]
    r_Earth2Sun = [1e6,  1.50147817e8,  2e6]
    assert sensors.isEclipse(r_sat, r_Earth2Sun, Re) == True
    
    r_sat = [600,7800,-500]
    r_Earth2Sun = [1e6,  1.50147817e8,  2e6]
    assert sensors.isEclipse(r_sat, r_Earth2Sun, Re) == False
    
    r_sat = [7000,-7800,-500]
    r_Earth2Sun = [1e7,  1.50147817e8,  2e8]
    assert sensors.isEclipse(r_sat, r_Earth2Sun, Re) == False
    
    r_sat = [7000,-7800,-500]
    r_Earth2Sun = [-1e7,  1.50147817e8,  1.4e8]
    assert sensors.isEclipse(r_sat, r_Earth2Sun, Re) == False
    
    r_sat = [6000,-6500,-6500]
    r_Earth2Sun = [1e7,  1.50147817e8,  1.4e8]
    assert sensors.isEclipse(r_sat, r_Earth2Sun, Re) == False

    r_sat = [5900,-6500,-6500]
    r_Earth2Sun = [1e7,  1.50147817e8,  1.4e8]
    assert sensors.isEclipse(r_sat, r_Earth2Sun, Re) == True