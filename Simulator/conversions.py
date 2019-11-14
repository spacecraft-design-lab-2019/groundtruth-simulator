# -*- coding: utf-8 -*-
import numpy as np
import math

#--------------------Coordinate Frames--------------------------
# TODO: Add conversion to/from NED (for use with pyIGRF)


def ECI_to_ECEF(r_ECI, GMST):
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
    lon = np.arctan2(r_ECEF[1], r_ECEF[0])
    alt = np.linalg.norm(r_ECEF) - rad_Earth

    return lat, lon, alt

#--------------------Miscellaneous-----------------------------

def skew(v):
    """
    Function: skew
        Calculates the skew matrix for cross-product calculation

    Inputs:

    """
    S = np.array([[0, -v[2], v[1]],
                 [v[2], 0, -v[0]],
                 [-v[1], v[0], 0]])
    return S


def mjd_2_GMST(mjd):
    # Reference: AA 279A Lecture 6, Slide 3
    d = mjd - 51544.5
    return math.fmod(np.radians(280.4606 + 360.9856473*d), 2*np.pi)

#--------------------Quaternions-----------------------------

def L(q):
    s = q[0]
    v = q[1:4]
    I3 = np.eye(3)

    upper = np.hstack([s, -v])
    lower = np.column_stack([v, s*I3 + skew(v)])
    return  np.row_stack([upper, lower])


def R(q):
    s = q[0]
    v = q[1:4]
    I3 = np.eye(3)

    upper = np.hstack([s, -v])
    lower = np.column_stack([v, s*I3 - skew(v)])
    return  np.row_stack([upper, lower])


def quat(v):
    assert len(v) == 4, f"can't make a quaternion out of: {v}"
    N = np.linalg.norm(v)
    assert N > 0, "Norm of input is 0. Can't make a quaternion out of this! Got {}".format(v)
    return v / N


def conj(q):
    q = np.copy(q)
    q[0] = -q[0]
    return q


def quatrot(q1, q2):
    vector = False
    if len(q2) == 3:
        q2 = np.append(0, q2)
        vector = True

    rotated = L(q1).T @ R(q1) @ q2

    if vector:
        return rotated[1:4]
    else:
        return rotated
