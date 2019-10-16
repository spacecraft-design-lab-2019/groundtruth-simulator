# -*- coding: utf-8 -*-
import numpy as np
import msise00
from datetime import datetime
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


def density_lookup(year,month,day,hour,altitude,glat,glon):
    """
    Function: density_lookup

    Gets atmospheric density using MSISE-00 atmospheric model.
    Must have https://pypi.org/project/msise00/ library installed.

    Inputs:
        year
        month
        day
        hour
        altitude (km)
        glat: geodetic latitude
        glon: geodetic longitude

    Outputs:
        rho: atmospheric density (kg/m^3)
    """
    atmos = msise00.run(time=datetime(year, month, day, hour), altkm=altitude, glat=glat, glon=glon)
    rho = atmos.Total.values[0].item()
    return rho

def dragCalc(r,v,cD,A,Re,wEarth,cmx,cmz,cpx,cpz,year,month,day,hour,altitude,glat,glon):
    R = np.linalg.norm(r)
    h = R - Re


    #constants for calculating density
    rho = density_lookup(year,month,day,hour,altitude,glat,glon)

    vRel = np.cross(wEarth,r)
    adrag = -0.5*rho*cD*A*norm(vRel)^2 * vRel/norm(vRel)

    #cp is center of pressure coordinate, cm is center of mass coorinate
    #note cp and cm must be in Local Vertical/Local Horizontal Coords
    mdrag = 0.5*rho*cD*A*norm(vRel)^2*np.array([cpx - cmx, 0, cpz - cmz])
    return adrag, mdrag


#-------------------------Torques--------------------------------

def gravityGradientTorque(r_sat, R_eci2principal, I, GM):
    """
    Function: gravityGradientTorque

    Inputs:
        r_sat: position vector of satellite in ECI
        R_eci2principal: rotation matrix from eci to principal axes
        I: moment of inertia matrix (in principal axes) (diagonal)
        GM: gravitational parameter of attracting body
    Outputs:
        M: moment due to gravity gradient
    """

    R = np.linalg.norm(r_sat)
    c = R_eci2principal @ r_sat
    
    Ix = I[0,0]
    Iy = I[1,1]
    Iz = I[2,2]
    cx = c[0]
    cy = c[1]
    cz = c[2]
    
    M = np.zeros((3,))
    M[0] = 3*GM/R**3 * (Iz - Iy)*cy*cx
    M[1] = 3*GM/R**3 * (Ix - Iz)*cz*cx
    M[2] = 3*GM/R**3 * (Iy - Ix)*cx*cy
    
    return M



#-------------------------------------------------------------
#-------------------------BACKUP STUFF -----------------------
#-------------------------------------------------------------

def accelHarmonic(state, R, n, m):
    # -------- TO FIX: BELOW CODE IS 1-INDEXED ------------

    x = 0
    y = 0
    z = 0
    V = np.zeros(50)
    W = np.zeros(50)

    r_vec = state[1:3]
    x_pos = r_vec[1]
    y_pos = r_vec[2]
    z_pos = r_vec[3]
    r = np.linalg.norm(r_vec)

    for j in range(0,m):
        for i in range(j,n):
            if i == 0:
                V[0,0] = R/r
                W[0,0] = 0
            elif i == j:
                V[i+1,j+1] = (2*j - 1)*((x_pos*R/r^2)*V[j,j] - (y_pos*R/r^2)*W[j,j])
                W[i+1,j+1] = (2*j - 1)*((x_pos*R/r^2)*W[j,j] - (y_pos*R/r^2)*V[j,j])
            elif i == j + 1:
                V[i+1,j+1] = ((2*i-1)/(i-j))*(z_pos*R/r^2)*V[i,j+1]
                W[i+1,j+1] = ((2*i-1)/(i-j))*(z_pos*R/r^2)*W[i,j+1]

            else:
                V[i+1,j+1] = ((2*i-1)/(i-j))*(z_pos*R/r^2)*V[i,j+1] - ((i+j-1)/(i-j))*(R^2/r^2)*V[i-1,j+1]
                W[i+1,j+1] = ((2*i-1)/(i-j))*(z_pos*R/r^2)*W[i,j+1] - ((i+j-1)/(i-j))*(R^2/r^2)*W[i-1,j+1]


    for j in range(0,m):
        for i in range(j,n):
            if j == 0:
                x = x + (GM/R^2)*(-C[i+1,1]*V[i+2,2]);
                y = y + (GM/R^2)*(-C[i+1,1]*W[i+2,2]);
            else:
                x = x + (GM/R^2)*0.5*((-C[i+1,j+1]*V[i+2,j+2] - S[i+1,j+1]*W[i+2,j+2]) + (np.math.factorial(i-j+2)/np.math.factorial(i-j))*(C[i+1,j+1]*V[i+2,j] + S[i+1,j+1]*W[i+2,j]))
                y = y + (GM/R^2)*0.5*((-C[i+1,j+1]*W[i+2,j+2] - S[i+1,j+1]*V[i+2,j+2]) + (np.math.factorial(i-j+2)/np.math.factorial(i-j))*(-C[i+1,j+1]*W[i+2,j] + S[i+1,j+1]*V[i+2,j]))

            z = z + (GM/R^2)*((i-j+1)*(-C[i+1,j+1]*V[i+2,j+1] - S[i+1,j+1]*W[i+2,j+1]))

    return np.array([x,y,z])

"""
function f_sphergrav = accelHarmonic2(state, GM, R, n, m)

X = load('constants.mat');
C = X.C;
S = X.S;

x = 0;
y = 0;
z = 0;
V = zeros(50);
W = zeros(50);

r_vec = state(1:3);
x_pos = r_vec(1);
y_pos = r_vec(2);
z_pos = r_vec(3);
r = norm(r_vec);

for j = 0:m
    for i = j:n
        if i == 0
            V(1,1) = R/r;
            W(1,1) = 0;

        elseif i == j
            V(i+1,j+1) = (2*j - 1)*((x_pos*R/r^2)*V(j,j) - (y_pos*R/r^2)*W(j,j));
            W(i+1,j+1) = (2*j - 1)*((x_pos*R/r^2)*W(j,j) - (y_pos*R/r^2)*V(j,j));

        elseif i == j + 1
            V(i+1,j+1) = ((2*i-1)/(i-j))*(z_pos*R/r^2)*V(i,j+1);
            W(i+1,j+1) = ((2*i-1)/(i-j))*(z_pos*R/r^2)*W(i,j+1);

        else
            V(i+1,j+1) = ((2*i-1)/(i-j))*(z_pos*R/r^2)*V(i,j+1) - ((i+j-1)/(i-j))*(R^2/r^2)*V(i-1,j+1);
            W(i+1,j+1) = ((2*i-1)/(i-j))*(z_pos*R/r^2)*W(i,j+1) - ((i+j-1)/(i-j))*(R^2/r^2)*W(i-1,j+1);
        end
    end
end

for j = 0:m
    for i = j:n
        if j == 0
            x = x + (GM/R^2)*(-C(i+1,1)*V(i+2,2));
            y = y + (GM/R^2)*(-C(i+1,1)*W(i+2,2));
        else
            x = x + (GM/R^2)*0.5*((-C(i+1,j+1)*V(i+2,j+2) - S(i+1,j+1)*W(i+2,j+2)) +...
                (factorial(i-j+2)/factorial(i-j))*(C(i+1,j+1)*V(i+2,j) + S(i+1,j+1)*W(i+2,j)));
            y = y + (GM/R^2)*0.5*((-C(i+1,j+1)*W(i+2,j+2) - S(i+1,j+1)*V(i+2,j+2)) +...
                (factorial(i-j+2)/factorial(i-j))*(-C(i+1,j+1)*W(i+2,j) + S(i+1,j+1)*V(i+2,j)));
        end

        z = z + (GM/R^2)*((i-j+1)*(-C(i+1,j+1)*V(i+2,j+1) - S(i+1,j+1)*W(i+2,j+1)));
    end
end

f_sphergrav = [x; y; z];
end

"""
