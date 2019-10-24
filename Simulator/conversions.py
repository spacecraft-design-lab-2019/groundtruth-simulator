# -*- coding: utf-8 -*-
import numpy as np

#--------------------Conversions-----------------------------

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
    long = np.arctan2(r_ECEF[1], r_ECEF[0])
    alt = np.linalg.norm(r_ECEF - rad_Earth)

    return lat, long, alt


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
    N = np.linalg.norm(v)
    assert N > 0, "Norm of v is 0. Can't make a quaternion out of this! Got {}".format(v)
    if len(v) == 3:
        return np.append(0, v) / N
    elif len(v) == 4:
        return v / N
    else:
        raise(Exception("can't make a quaternion out of: {}".format(v)))

def conj(q):
    q = copy(q)
    q[0] = -q[0]
    return q

def quatrot(q1, q2):
    if len(q2) == 3:
        vector = True
    q1 = quat(q1)
    q2 = quat(q2)

    rotated = L(q1).T @ R(q1) @ q2

    if vector:
        return rotated[1:4]
    else:
        return rotated
