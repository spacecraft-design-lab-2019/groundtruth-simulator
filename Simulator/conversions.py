# -*- coding: utf-8 -*-
import numpy as np
import math

#--------------------Coordinate Frames--------------------------

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


def ECEF_to_ECI(r_ECEF, GMST):
    """
    Function: ECEF_to_ECI
        Converts position vector in ECEF to ECI.

    Inputs:
        r_ECEF: position vector in Earth Centered Earth Fixed (ECEF)
        GMST: current Greenwich Mean Sidereal Time [rad]

    Outputs:
        r_ECI: position vector in Earth Centered Inertial (ECI)
    """

    rotation = np.array([[np.cos(GMST), np.sin(GMST), 0],
                         [-np.sin(GMST), np.cos(GMST), 0],
                         [0, 0, 1]])
    r_ECI = rotation.transpose() @ r_ECEF

    return r_ECI


def ECEF_to_LLA(r_ECEF, rad_Earth):
    """
    Function: ECEF_to_LLA
        Converts position vector in ECEF to geocentric coordinates.

    Inputs:
        r_ECEF: position vector in Earth Centered Earth Fixed (ECEF)
        rad_Earth: radius of the Earth [km]

    Outputs:
        lat:    latitude    [rad]
        long:   longitude   [rad]
        alt:    altitude    [rad]
    """

    # lat = np.arcsin(r_ECEF[2] / np.linalg.norm(r_ECEF))
    lat = np.arcsin(r_ECEF[2] / norm2(r_ECEF))
    lon = np.arctan2(r_ECEF[1], r_ECEF[0])
    # alt = np.linalg.norm(r_ECEF) - rad_Earth
    alt = norm2(r_ECEF) - rad_Earth

    return np.array([lat, lon, alt])


def NED_to_ECI(vec_NED, glat, glon, GMST):
    """
    Function: NED_to_ECI
        Converts a vector from NED to ECI coordinates.

    Inputs:
        vec_NED: position vector in NED coordinates
        glat: geodetic latitude (radians)
        glon: geodetic longitude (radians)
        GMST: Greenwich mean sidereal time (radians)
    """
    R = np.array([[-np.sin(glat)*np.cos(glon), -np.sin(glon), -np.cos(glat)*np.cos(glon)],
                  [-np.sin(glat)*np.sin(glon),  np.cos(glon), -np.cos(glat)*np.sin(glon)],
                  [              np.cos(glat),           0.0,              -np.sin(glat)]])

    vec_ECEF = R @ vec_NED
    vec_ECI = ECEF_to_ECI(vec_ECEF, GMST)
    return vec_ECI



#--------------------Quaternions-----------------------------

def L(q):
    """
    Matrix equivilent to left-hand quaternion for mulitplication
    (eg. q1 * q2 = L(q1)q2)
    """
    s = q[0]
    v = q[1:4]
    L = np.array([[s, -v[0], -v[1], -v[2]],
                  [v[0], s, -v[2], v[1]],
                  [v[1], v[2], s, -v[0]],
                  [v[2], -v[1], v[0], s]])
    return L


def R(q):
    """
    Matrix equivilent to right-hand quaternion for mulitplication
    (eg. q1 * q2 = R(q2)q1)
    """
    s = q[0]
    v = q[1:4]
    R = np.array([[s, -v[0], -v[1], -v[2]],
                  [v[0], s, v[2], -v[1]],
                  [v[1], -v[2], s, v[0]],
                  [v[2], v[1], -v[0], s]])
    return R


def quat(v):
    """
    Normalizes a quaternion (or any 4d vector)
    """
    assert len(v) == 4, f"can't make a quaternion out of: {v}"
    N = np.linalg.norm(v)
    assert N > 0, "Norm of input is 0. Can't make a quaternion out of this! Got {}".format(v)
    return v / N

def conj(q):
    """
    Calculates the quaternion conjugate
    """
    q = np.copy(q)
    q[0] = -q[0]
    return q

def quatrot(q1, q2):
    """
    Rotates a vector `q2` using the quaternion `q1`.
    """
    vector = False
    if len(q2) == 3:
        q2 = np.append(0, q2)
        vector = True

    rotated = L(q1) @ R(q1).T @ q2

    if vector:
        return rotated[1:4]
    else:
        return rotated

# def quatmult(q1, q2):
#     """
#     Multiplies two quaternions as in q1*q2 -> L(q1)*q2.
#     """
#     return L(q1) @ q2


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
    """
    Function: mjd_2_GMST
        Calculates Greenwich Mean Sidereal Time

    Inputs:
        mjd - modified julian day
    Outputs:
        GMST - Greenwich Mean Sidereal Time [rad]

    Reference: Vallado
    """
    d = (mjd - 51544.5) / 36525.0
    GMST = 67310.54841 + (876600*3600 + 8640184.812866) * d + 0.093104 * d**2 - 6.2 * 10**(-6)*d**3
    return (GMST % 86400) / 240 * (np.pi / 180)
    # d = mjd - 51544.5
    # return math.fmod(np.radians(280.4606 + 360.9856473*d), 2*np.pi)

def unit(ax):
    """
    make a unit vector from int or str input (0, 1, 2) / ('x', 'y', 'z')
    """
    if isinstance(ax, str):
        if ax == 'x':
            ax = 0
        elif ax == 'y':
            ax = 1
        elif ax == 'z':
            ax = 2
        else:
            raise(Exception("invalid unit axis input {}".format(ax)))

    N = np.zeros(3)
    N[ax] = 1
    return N

def cross3(v1, v2):
    """
    Calculates the cross-product of a pair of vectors in R^3
    """
    x = ((v1[1] * v2[2]) - (v1[2] * v2[1]))
    y = ((v1[2] * v2[0]) - (v1[0] * v2[2]))
    z = ((v1[0] * v2[1]) - (v1[1] * v2[0]))

    return np.array([x, y, z])


def norm2(v):
    """
    Caluclates the 2-norm of a vector
    """
    # return (v.T @ v) ** (0.5)
    return math.sqrt(sum(x*x for x in v))