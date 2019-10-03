# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 14:24:20 2019

@author: Kevin
"""
import numpy as np
#from numpy import linalg as LA
def dragCalc(r,v,M,A,h0,rho0,H,Re,wEarth):
    R = np.linalg.norm(r)
    h = R - Re
    vRel = np.cross(wEarth,r)
    rho = rho0 * exp(-(h-h0)/H)
    B = Cd*A/M
    adrag = -0.5*B*rho*norm(vRel)^2 * vRel/norm(vRel)
    return adrag
    
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