# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 17:22:15 2019

"""
import pytest
import os, sys, inspect
import sun_sensor_math as LA

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)

def pytest_addoption(parser):
    parser.addoption(
        "--cmdopt", action="store", default="type1", help="my option: type1 or type2"
    )

@pytest.fixture
def cmdopt(request):
    return request.config.getoption("--cmdopt")

def test_sense2vector():
    assert True
    
def test_vector2sense():
    assert True
    
def test_isEclipse2():
    measurements = [3,3,3,3,3,3]
    thresh = 5
    assert LA.isEclipse2(measurements, thresh) == True
    
    measurements = [6,5,6,5,6,7]
    thresh = 5
    assert LA.isEclipse2(measurements, thresh) == False
    
    measurements = [6.1,101,198,32,.1,7]
    thresh = 5
    assert LA.isEclipse2(measurements, thresh) == False
    
    measurements = [6.14,1861,198,32.5,.1,199]
    thresh = 200.1
    assert LA.isEclipse2(measurements, thresh) == True
    
    measurements = [6.14,1861,198,32.5,.1,200.1]
    thresh = 200.1
    assert LA.isEclipse2(measurements, thresh) == False