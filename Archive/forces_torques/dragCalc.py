# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 14:24:20 2019

"""
import numpy as np
import msise00
from datetime import datetime
#from numpy import linalg as LA

#atmospheric density
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

#testing
print(density_lookup(2019,10,6,12,500,90,90))

def dragCalc(r,v,cD,A,Re,wEarth,cmx,cmz,cpx,cpz,year,month,day,hour,altitude,glat,glon):
    #constants for calculating density
    rho = density_lookup(year,month,day,hour,altitude,glat,glon)

    vRel = np.cross(wEarth,r)
    adrag = -0.5*rho*cD*A*np.linalg.norm(vRel)^2 * vRel/np.linalg.norm(vRel)

    #cp is center of pressure coordinate, cm is center of mass coorinate
    #note cp and cm must be in Local Vertical/Local Horizontal Coords
    mdrag = 0.5*rho*cD*A*np.linalg.norm(vRel)^2*np.array([cpx - cmx, 0, cpz - cmz])
    return adrag, mdrag

"""
function adrag = dragCalc(r,v)

R = norm(r);
M = 1500; %kg
Cd = 2.3;
A = 20e-6; %m^2
h0 = 0;
rho0 = 1.225e9;
H = 10; %km
Re = 6378;
h = R - Re;
wEarth = [0 0 2*pi/86164.1];
vRel = v - cross(wEarth,r);

rho = rho0 * exp(-(h-h0)/H);
B = Cd*A/M;
fdrag = -0.5*B*rho*norm(vRel)^2 * vRel/norm(vRel);
adrag = fdrag;

end
"""
