# -*- coding: utf-8 -*-
import numpy as np

#---------------------------Dynamics---------------------------------

def accelPointMass(r_sat, r_body, GM):
    #--------currently Earth gravity only assuming r_sat is in ECI
    #----------will eventually change to any body given r_body
    accel = -GM * r_sat / (np.linalg.norm(r_sat)**3)
    return accel


def accelHarmonic():
    return 0


