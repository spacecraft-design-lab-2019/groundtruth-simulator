# -*- coding: utf-8 -*-
import numpy as np
import constants as con


#-------------------------Forces---------------------------------

def gravityPointMass(r_sat, r_body, GM):
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


def gravityEarthJ2(r_sat, GM, J2, rad_Earth):
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
    
    f = np.zeros((3,))
    f = -0.5*GM*J2*rad_Earth**2 * (3/r**4 - 15*(z**2 / r**6)) * R
    f[2] = f[2] + -0.5*GM*J2*rad_Earth**2 * (6 * (z / r**5))
    
    return f


def dragCalc(r_ECI,v_ECI,GMST, mjd, cD,A,cmx,cmz,cpx,cpz):
    #constants for calculating density
    
    rho = con.Environment.density_lookup(r_ECI, GMST, mjd, con.Earth.radius)
    vRel = np.cross(con.Earth.w, r_ECI)
    adrag = -0.5*rho*cD*A*np.linalg.norm(vRel)^2 * vRel/np.linalg.norm(vRel)

    #cp is center of pressure coordinate, cm is center of mass coorinate
    #note cp and cm must be in Local Vertical/Local Horizontal Coords
    A = con.SpacecraftStructure.surfArea
    cD = con.SpacecraftStructure.cD
    mdrag = 0.5*rho*cD*A*np.linalg.norm(vRel)^2*np.array([cpx - cmx, 0, cpz - cmz])
    return adrag, mdrag


#-------------------------Torques--------------------------------

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

