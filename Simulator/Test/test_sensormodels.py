import pytest
import os, sys, inspect

# add current folder to the path
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)

from sensormodels import *
import sim_config as config
import numpy as np

def test_init():

    # constructors shouldn't error obviously:
    S = SpacecraftSensors(config.mag_params, config.gyro_params, config.sun_params)

    x = np.random.rand(3)

    # and neither should these:
    S.gyroscope.measure(x)
    S.magnetometer.measure(x)
    S.sunsensor.measure(x)

    # and this too...
    S.gyroscope.update()
    S.gyroscope.measure(x)

    # identity sensor:
    lme = LinearErrorModel.withDim(3)
    np.testing.assert_allclose(lme.measure(x), x,  atol = 1e-10)

    # bias only sensor
    lme.b = 1
    np.testing.assert_allclose(lme.measure(x), x+1,  atol = 1e-10)

    lme.b = np.ones(3)
    np.testing.assert_allclose(lme.measure(x), x+1,  atol = 1e-10)

    lme.T = np.eye(3)
    np.testing.assert_allclose(lme.measure(x), (2*x)+1,  atol = 1e-10)