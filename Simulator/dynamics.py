# -*- coding: utf-8 -*-
import numpy as np
from conversions import *
from numpy import linalg as LA
from constants import *
#-------------------------Dynamics---------------------------------

def gravityPointMass(r_sat, GM, r_body=np.zeros((3,))):
    """
    Function: gravityPointMass
        Calculates the gravitational force of a rigid spherical body.

    Inputs:
        r_sat: position vector of the satellite in ECI
        r_body: position vector of the attracting body in ECI
        GM: gravitaional parameter of attracting body

    Outputs:
        accel: acceleration due to gravity
    """

    accel = -GM * (r_sat - r_body) / (np.linalg.norm(r_sat - r_body)**3)
    return accel


def gravityEarthJ2(r_sat, earth):
    """
    Function: gravityEarthJ2

    Inputs:
        r_sat: position vector of satellite in ECI
        GM: gravitational parameter of Earth
        J2: constant representing non-spherical shape of Earth
        rad_Earth: radius of the Earth
    Outputs:
        f: acceleration due to Earth

    Reference: AA 279A, Lecture 14, Slide 5
    """

    z = r_sat[2]
    r = np.linalg.norm(r_sat)
    R = r_sat / r

    f = np.zeros(3)
    f = -0.5*earth.GM*earth.J2*earth.radius**2 * (3/r**4 - 15*(z**2 / r**6)) * R
    f[2] = f[2] + -0.5*earth.GM*earth.J2*earth.radius**2 * (6 * (z / r**5))

    return f


def dragCalc(state, environment, structure):

    # in ECI:
    r = state[0:3]
    v = state[7:10]
    # q is body to ECI
    q = state[3:7]

#    rho = environment.density_lookup(r, GMST(mjd), mjd)
    #drag equation fit coefficients
    a = 4.436e-09
    b = -0.01895
    c = 4.895e-12
    d = -0.008471
    R = LA.norm(r) - environment.earth.radius
    rho = a*np.exp(b*R) + c*np.exp(d*R)

    vRel = v - np.cross(environment.earth.w, r)
    vRel_body = quatrot(conj(q), vRel) # get from inertial to body

    return structure.aerodrag(rho, vRel_body)


def gravityGradientTorque(r_sat, I, GM):
    """
    Function: gravityGradientTorque

    Inputs:
        r_sat: position vector of satellite in ECI
        I: moment of inertia matrix
        GM: gravitational parameter of attracting body
    Outputs:
        M: moment due to gravity gradient
    """

    return (3*GM / ((r_sat.T @ r_sat)**5/2)) * np.cross(r_sat, I @ r_sat)

