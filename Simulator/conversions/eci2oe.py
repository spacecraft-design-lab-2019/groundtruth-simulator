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


    eVec = (1/GM)*((v^2 - GM/r)*r_ijk - dot(r_ijk, v_ijk)*v_ijk);
    e = norm(eVec);

    #Compute the size of the orbit
    mechEnergy = 0.5*v^2 - GM/r

    if e !=1:
        a = -GM/(2*mechEnergy)
        p = a*(1 - e^2)
    else:
        a = inf
        p = h^2/GM

# Compute the orientation of the orbit
i = acosd(hVec(3)/h);
Omega = np.arccos(nVec[0]/n);
omega = np.arccos(np.dot(nVec, eVec)/(n*e));
nu = acosd(dot(eVec, r_ijk)/(e*r));

# Place angles in the correct domains
if nVec(2) < 0:
    Omega = 360 - Omega

if eVec(3) < 0:
    omega = 360 - omega

if dot(r_ijk, v_ijk) < 0:
    nu = 360 - nu

# Account for any special cases
if i == 0 and e != 0:  #Elliptical equatorial
    # Provide the longitude of periapsis (PI = Om + w)
    ang = acosd(eVec(1)/e);
    
    if eVec(2) < 0
        ang = 360 - ang;
    end
elseif i~= 0 && e == 0 % Circular inclined
    % Provide the argument of latitude (u = w + anom)
    ang = acosd(dot(nVec, r_ijk)/(n*r));
    
    if r_ijk(3) < 0
        ang = 360 - ang;
    end
elseif i == 0 && e == 0 % Circular equatorial
    % Provide the true latitude (lambda = Om + w + anom)
    ang = acosd(r_ijk(1)/r);
    
    if r_ijk(2) < 0
        ang = 360 - ang;
    end
else
    % Default output for ang
    ang = NaN;
end

oe(1) = a;
oe(2) = e;
oe(3) = i;
oe(4) = Omega;
oe(5) = omega;
oe(6) = nu;
oe(7) = ang;

oe(3:7) = deg2rad(oe(3:7));
end

"""
