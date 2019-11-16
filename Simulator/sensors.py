# -*- coding: utf-8 -*-
"""
Integration of Sensor Model
"""

import numpy as np
from sensormodels import *
from sun_model import *


class SpacecraftSensors():
    """
    A class to store spacecraft sensors
    """
    def __init__(self):
        self.magnetometer = Sensor(errormodel=LinearErrorModel(scaleF=mag_params["scaleF"], caSense=mag_params["caSense"], b=mag_params["b"], cov=mag_params["cov"]))
        self.gyroscope    = Sensor(errormodel=LinearErrorModel(scaleF=gyro_params["scaleF"], caSense=gyro_params["caSense"], b=gyro_params["b"], cov=gyro_params["cov"]))
        self.sunsensor    = Sensor(errormodel=LinearErrorModel(scaleF=sun_params["scaleF"], caSense=sun_params["caSense"], b=sun_params["b"], cov=sun_params["cov"]))


def sunVector(rsat,jdate):
    """
    Inputs:
        rsat: ECI position of satellite
        jdate: mean julian date time
    Output:
        position vector (3-vector) from satellite to sun
    """
    earth2sun = sun1(jdate)
    rsun = (earth2sun - rsat)/LA.norm(earth2sun - rsat)
    return rsun