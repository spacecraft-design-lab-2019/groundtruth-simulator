from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_sun
import astropy
import numpy as np
import math

def approx_sun_position_ECI(MJD):
    """
    This is using the equations given in Motenbruck and Gill's Satellite Orbits book
    Inputs:
    MJD - Modified Julian Day (J2000) as a Real Number
    Outputs:
    r_vec - numpy array with x, y, z of Sun position in ECI at input time
    """
    JD = MJD + 2400000.5
    OplusW = 282.94
    T = (JD - 2451545.0) / 36525

    M = math.radians(357.5256 + 35999.049 * T)

    lon = math.radians(OplusW + math.degrees(M) + 6892 / 3600 * math.sin(M) + 72 / 3600 * math.sin(2*M))
    r_mag = (149.619 - 2.499 * math.cos(M) - 0.021 * math.cos(2*M)) * 10**6

    epsilon = math.radians(23.43929111)
    r_vec = (r_mag * math.cos(lon), r_mag * math.sin(lon) * math.cos(epsilon), r_mag * math.sin(lon) * math.sin(epsilon))

    return r_vec

def sun_position_ECI(dt):
    """
    Queries astropy module for ECI sun position.

    Inputs:
        dt - datetime object
    Outputs:
        rECI - position vector to the sun in ECI coordinates (km)
    """
    AU_km = 149597870.7
    t = Time(dt, scale='utc', format='datetime')
    sun = get_sun(t)

    r_GCRS = np.array([sun.ra.value, sun.dec.value, sun.distance.value * AU_km])
    r_ECI = GCRS_to_ECI(r_GCRS)
    return r_ECI


def GCRS_to_ECI(r_GCRS):
    """
    Converts from GCRS to ECI.

    Inputs:
        r_GCRS - (ra, dec, distance) in (deg, deg, km)
    Outputs:
        r_ECI - position vector in ECI (km)
    """

    ra = r_GCRS[0] * np.pi / 180 # (rad)
    dec = r_GCRS[1] * np.pi / 180 # (rad)
    dist = r_GCRS[2] # km

    r_ECI = np.array([dist * np.cos(ra),
                      dist * np.sin(ra) * np.cos(-dec),
                      dist * np.sin(ra) * np.sin(-dec)])

    #theta = (280.4606 + 360.9856473 * (MJD - 51544.5))/180 * np.pi
    #Rz = np.array([[np.cos(theta), np.sin(theta),0],
    #    [-np.sin(theta), np.cos(theta), 0],
    #    [0,0,1]])

    return r_ECI



# MJD = 53005
# t = Time(MJD, scale='utc', format='mjd')
# AU_km = 149597870.7
# sun = get_sun(t)
# r_GCRS = np.array([sun.ra.value, sun.dec.value, sun.distance.value * AU_km])
# r_ECI = GCRS_to_ECI(r_GCRS,MJD)
