# -*- coding: utf-8 -*-
"""
Function: eci2oe

Converts ECI position and velocity to orbital elements.
"""
#test code
import numpy as np
from numpy import linalg as LA
state = np.array([1,2,3,4,5,6])

def eci2oe(state, GM):
#state is a 6x1 vector of position and velocity    
    print("hello")
    r_ijk = state[0:3]
    v_ijk = state[3:6]
    r = LA.norm(r_ijk)
    v = LA.norm(v_ijk)

    #Create all necessary vectors
    hVec = np.cross(r_ijk, v_ijk)
    h = LA.norm(hVec)

    nVec = np.cross(np.array([0, 0, 1]), hVec)
    n = LA.norm(nVec)


    eVec = (1/GM)*((v^2 - GM/r)*r_ijk - np.dot(r_ijk, v_ijk)*v_ijk)
    e = LA.norm(eVec)

    #Compute the size of the orbit
    mechEnergy = 0.5*v^2 - GM/r

    if e !=1:
        a = -GM/(2*mechEnergy)
        p = a*(1 - e^2)
    else:
        a = np.inf
        p = h^2/GM

    # Compute the orientation of the orbit
    i = np.arccos(hVec(3)/h)
    Omega = np.arccos(nVec[0]/n)
    omega = np.arccos(np.dot(nVec, eVec)/(n*e))
    nu = np.arccos(np.dot(eVec, r_ijk)/(e*r))
    
    # Place angles in the correct domains
    if nVec(2) < 0:
        Omega = 2*np.pi - Omega
    
    if eVec(3) < 0:
        omega = 2*np.pi - omega
    
    if np.dot(r_ijk, v_ijk) < 0:
        nu = 2*np.pi - nu
    
    # Account for any special cases
    if i == 0 and e != 0:  #Elliptical equatorial
        # Provide the longitude of periapsis (PI = Om + w)
        ang = np.arccos(eVec(1)/e)
        
        if eVec(2) < 0:
            ang = 2*np.pi - ang
    elif i != 0 and e == 0: # Circular inclined
        # Provide the argument of latitude (u = w + anom)
        ang = np.arccos(np.dot(nVec, r_ijk)/(n*r))
        
        if r_ijk(3) < 0:
            ang = 2*np.pi - ang
            
    elif i == 0 and e == 0: # Circular equatorial
        # Provide the true latitude (lambda = Om + w + anom)
        ang = np.arccos(r_ijk[0]/r)
        
        if r_ijk(2) < 0:
            ang = 2*np.pi - ang
    else:
        # Default output for ang
        ang = np.nan
        
    oe = np.zeros((6,1))
    oe[0] = a
    oe[1] = e
    oe[2] = i
    oe[3] = Omega
    oe[4] = omega
    oe[5] = nu
    oe[6] = ang


