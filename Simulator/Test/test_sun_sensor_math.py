# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 17:22:15 2019

"""
import pytest
import os, sys, inspect
import sun_sensor_math as LA
import numpy as np

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)

def test_sense2vector_error():
    with pytest.raises(TypeError):
        measurements = 5
        r_Earth2Sun = [1,1,1]
        r_sat = [1,1,1]
        assert LA.sense2vector(measurements, r_Earth2Sun, r_sat)

    with pytest.raises(TypeError):        
        measurements = [5,3,3,4,5,6]
        r_Earth2Sun = 1
        r_sat = [1,1,1]
        assert LA.sense2vector(measurements, r_Earth2Sun, r_sat)

    with pytest.raises(TypeError):        
        measurements = [5,3,3,4,5,6]
        r_Earth2Sun = np.array([1,1,1])
        r_sat = [1,1,1]
        assert LA.sense2vector(measurements, r_Earth2Sun, r_sat)

    with pytest.raises(TypeError):        
        measurements = [5,3,3,4,5,6]
        r_Earth2Sun = [1,1,1]
        r_sat = 1
        assert LA.sense2vector(measurements, r_Earth2Sun, r_sat)
    
def test_vector2sense():
    assert True
    
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

    