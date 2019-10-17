# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import os

#directory
dir_path = os.path.dirname(os.path.realpath(__file__))
mSun = 1.98855e30 #kg
GM = 3.986e5


class Earth():
    """
    A class to store Earth parameters
    """
    def __init__(self,
                R = 6378.1378366, #Earth Equatorial Radius
                w = 2*np.pi/86164.1, #orbital speed earth
                mass = 5.9742e24, #kg
                SMA = 149598023, #semimajor axis, km
                GM = 3.986e5): #gravitational param, km^3/s^2            
        self.R = R
        self.w = w
        self.mass = mass
        self.SMA = SMA
        self.GM = GM

class Moon():
    """
    A class to store Moon parameters
    """
    def __init__(self,
                R = 1738, #Equatorial Radius, km
                mass = 7.3483e22, #kg
                SMA = 38400, #semimajor axis, km
                GM = 4902.799): #gravitational param, km^3/s^2            
        self.R = R
        self.mass = mass
        self.SMA = SMA
        self.GM = GM