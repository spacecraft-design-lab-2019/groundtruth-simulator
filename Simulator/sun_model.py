from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_sun
import astropy
import numpy as np
import math
import julian

def approx_sun_position_ECI(mjd):
    """
    This is using the equations given in Motenbruck and Gill's Satellite Orbits book
    Inputs:
    mjd - Modified Julian Day (J2000) as a Real Number
    Outputs:
    r_vec - numpy array with x, y, z of Sun position in ECI at input time
    """
    JD = mjd + 2400000.5
    OplusW = 282.94
    T = (JD - 2451545.0) / 36525

    M = math.radians(357.5256 + 35999.049 * T)

    lon = math.radians(OplusW + math.degrees(M) + 6892 / 3600 * math.sin(M) + 72 / 3600 * math.sin(2*M))
    r_mag = (149.619 - 2.499 * math.cos(M) - 0.021 * math.cos(2*M)) * 10**6

    epsilon = math.radians(23.43929111)
    r_vec = (r_mag * math.cos(lon), r_mag * math.sin(lon) * math.cos(epsilon), r_mag * math.sin(lon) * math.sin(epsilon))

    return r_vec

def sun_position_ECI(mjd):
    """
    Queries astropy module for ECI sun position.

    Inputs:
        mjd - modified julian day
    Outputs:
        rECI - position vector to the sun in ECI coordinates (km)
    """
    AU_km = 149597870.7
    t = Time(mjd, format='mjd')
    sun = get_sun(t).cartesian
    
    r_ECI = AU_km * np.array([sun.x.value, sun.y.value, sun.z.value])
    return r_ECI
