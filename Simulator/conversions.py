# -*- coding: utf-8 -*-
import numpy as np

#--------------------Conversions-----------------------------

def ECI_to_ECEF(r_ECI, GMST)
    """
    Function: ECI_to_ECEF
        Converts position vector in ECI to ECEF.
        
    Inputs:
        r_ECI: position vector in Earth Centered Inertial (ECI)
        GMST: current Greenwich Mean Sidereal Time [rad]
        
    Outputs:
        r_ECEF: position vector in Earth Centered Earth Fixed (ECEF)
    """
    
    rotation = np.array([[np.cos(GMST), np.sin(GMST), 0],
                         [-np.sin(GMST), np.cos(GMST), 0],
                         [0, 0, 1]])
    r_ECEF = rotation @ r_ECI
    
    return r_ECEF


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

