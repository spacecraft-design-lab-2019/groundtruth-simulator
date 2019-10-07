# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 14:24:20 2019

"""
import numpy as np
#from numpy import linalg as LA
def dragCalc(r,v,cD,A,Re,wEarth,cmx,cmz,cpx,cpz):
    R = np.linalg.norm(r)
    h = R - Re
    
    
    #constants for calculating density
    a = 4.436e-09;
    b = -0.01895;
    c = 4.895e-12;
    d = -0.008471;
    rho = a*exp(b*h) + c*exp(d*h) #based on JD-2008 model
    
    vRel = np.cross(wEarth,r)
    adrag = -0.5*rho*cD*A*norm(vRel)^2 * vRel/norm(vRel)
    
    #cp is center of pressure coordinate, cm is center of mass coorinate
    #note cp and cm must be in Local Vertical/Local Horizontal Coords
    mdrag = 0.5*rho*cD*A*norm(vRel)^2*np.array([cpx - cmx, 0, cpz - cmz]) 
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