# -*- coding: utf-8 -*-
import pdb
import numpy as np
from conversions import *
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
    f[2] += -0.5*earth.GM*earth.J2*earth.radius**2 * (6 * (z / r**5))

    return f


def dragCalc(state, environment, structure):

    r = state[0:3]
    v = state[7:10]
    q = state[3:7]

    rho = environment.density_lookup(r)

    vRel = v - np.cross(environment.earth.w, r)
    vRel_body = quatrot(conj(q), vRel) # get from inertial to body

    # Simplified model - assumes constant wetted area:
    # A = 0.01
    # B = structure.cD*A/structure.mass
    # adrag = -0.5*B*rho*np.linalg.norm(vRel) * vRel
    # mdrag = np.zeros(3)
    # pdb.set_trace()
    # return adrag, mdrag

    adrag, mdrag = structure.aerodrag(rho, vRel_body)
    return quatrot(q, adrag), quatrot(q, mdrag)


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

