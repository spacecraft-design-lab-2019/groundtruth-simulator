# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 08:59:18 2016

Demo of particulate contaminant transport on an orbiting satellite with solar drag

https://www.particleincell.com/2016/particulate-orbit/

Equation for solar position based on code by David Eagle from:
http://www.mathworks.com/matlabcentral/fileexchange/39703-cowells-method-for-earth-satellites/content/sun1.m
"""

import numpy as np
import math
# import pylab as pl
# import sys

#start equations from computing solar position from code by David Eagle
def r2r (x):
    # revolutions to radians function       
    #  x = argument (revolutions; 0 <= x <= 1)    
    #  y = equivalent x (radians; 0 <= y <= 2 pi)    
    y = 2.0 * np.pi * (x - np.fix(x));
    return y
    
def atan3 (a, b):
    # four quadrant inverse tangent
    #  a = sine of angle
    #  b = cosine of angle
    #  y = angle (radians; 0 =< c <= 2 * pi)   
    epsilon = 0.0000000001;
    pidiv2 = 0.5 * np.pi;
    
    if (abs(a) < epsilon):
       y = (1 - np.sign(b)) * pidiv2;
       return y
    else:
       c = (2 - np.sign(a)) * pidiv2;
    
    if (abs(b) < epsilon):
       y = c;
       return y
    else:
       y = c + np.sign(a) * np.sign(b) * (abs(np.arctan(a / b)) - pidiv2);
    return y
    
def sun1(jdate):
# solar ephemeris
#  jdate = julian day
# output
#  rasc = right ascension of the sun (radians)         (0 <= rasc <= 2 pi)
#  decl = declination of the sun (radians)        (-pi/2 <= decl <= pi/2)
#  rsun = eci position vector of the sun (meters)

# note, coordinates are inertial, geocentric, equatorial and true-of-date

    atr = np.pi / 648000;
    
    rsun = np.zeros(3);
    
    # time arguments
    
    #CE 2000 January 01 12:00:00.0 UT  Saturday
    djd = jdate - 2451545;
    
    t = (djd / 36525) + 1;
    
    #fundamental arguments (radians)    
    gs = r2r(0.993126 + 0.0027377785 * djd);
    lm = r2r(0.606434 + 0.03660110129 * djd);
    ls = r2r(0.779072 + 0.00273790931 * djd);
    g2 = r2r(0.140023 + 0.00445036173 * djd);
    g4 = r2r(0.053856 + 0.00145561327 * djd);
    g5 = r2r(0.056531 + 0.00023080893 * djd);
    rm = r2r(0.347343 - 0.00014709391 * djd);
    
    # geocentric, ecliptic longitude of the sun (radians)   
    plon = 6910 * np.sin(gs) + 72 * np.sin(2 * gs) - 17 * t * np.sin(gs);
    plon = plon - 7 * np.cos(gs - g5) + 6 * np.sin(lm - ls) + 5 * np.sin(4 * gs - 8 * g4 + 3 * g5);
    plon = plon - 5 * np.cos(2 * (gs - g2)) - 4 * (np.sin(gs - g2) - np.cos(4 * gs - 8 * g4 + 3 * g5));
    plon = plon + 3 * (np.sin(2 * (gs - g2)) - np.sin(g5) - np.sin(2 * (gs - g5)));
    plon = ls + atr * (plon - 17 * np.sin(rm));
    
    # geocentric distance of the sun (km)   
    rsm = 149597870.691 * (1.00014 - 0.01675 * np.cos(gs) - 0.00014 * np.cos(2 * gs));
   
    #change to meters
    rsm *= 1000
    
    # obliquity of the ecliptic (radians)  
    obliq = atr * (84428 - 47 * t + 9 * np.cos(rm));
    
    # geocentric, equatorial right ascension and declination (radians)
    a = np.sin(plon) * np.cos(obliq);
    b = np.cos(plon);
       
#    rasc = atan3(a, b);
    rasc = math.atan2(a, b);
    decl = np.arcsin(np.sin(obliq) * np.sin(plon));
    
    # geocentric position vector of the sun (kilometers)
    rsun[0] = rsm * np.cos(rasc) * np.cos(decl);
    rsun[1] = rsm * np.sin(rasc) * np.cos(decl);
    rsun[2] = rsm * np.sin(decl);
    
    return rsun
#
#rsun = sun1(1929)
#print(rsun)
