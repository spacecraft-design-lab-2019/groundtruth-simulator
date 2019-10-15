# -*- coding: utf-8 -*-
import numpy as np

#--------------------Conversions-----------------------------

def ECEF_to_LLA(r_ECEF, rad_Earth):
    """
    Function: ECEF_to_LLA
        Converts position vector in ECEF to geodetic coordinates.
        
    Inputs:
        r_ECEF: position vector in Earth Centered Earth Fixed (ECEF)
        rad_Earth: radius of the Earth [km]
        
    Outputs:
        lat:    latitude    [rad]
        long:   longitude   [rad]
        alt:    altitude    [rad]
    """
    
    lat = np.arcsin(r_ECEF[2] / np.linalg.norm(r_ECEF))
    long = np.arctan2(r_ECEF[1], r_ECEF[0])
    alt = np.linalg.norm(r_ECEF - rad_Earth)
    
    return lat, long, alt