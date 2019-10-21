# -*- coding: utf-8 -*-
import numpy as np
# import msise00


#--------------------Spacecraft Structure----------------------------

class SpacecraftStructure():
    """
    A class to store spacecraft structural poperties
    """
    def __init__(self,
                 I = np.array([[17,0,0],[0,18,0],[0,0,22]]),
                 surfArea = np.array([.0025, 9.28e-5,1.024e-4,.003712]),
                 normVec1 = np.array([1,0,0]),
                 normVec2 = np.array([0,1,0]),
                 normVec3 = np.array([0,0,1]),
                 normVec4 = np.array([-1,0,0]),
                 normVec5 = np.array([0,-1,0]),
                 normVec6 = np.array([0,0,-1]),
                 cD = 2.3): #drag coefficient
        self.I = I
        self.surfArea = surfArea
        self.normVec1 = normVec1
        self.normVec2 = normVec2
        self.normVec3 = normVec3
        self.normVec4 = normVec4
        self.normVec5 = normVec5
        self.normVec6 = normVec6
        self.cD = cD


#-------------------------Environment---------------------------------

class Environment():
    """
    A class to store environment constants / lookup functions.
    """
    def __init__(self):
        self.earth = Earth()

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


class Earth():
    """
    A class to store Earth parameters
    """
    def __init__(self,
                radius = 6378.1378366, #Earth Equatorial Radius
                w = 2*np.pi/86164.1, #orbital speed earth
                mass = 5.9742e24, #kg
                SMA = 149598023, #semimajor axis, km
                J2 = 1.0826e-3, # J2 constant
                GM = 3.986e5): #gravitational param, km^3/s^2            
        self.radius = radius
        self.w = w
        self.mass = mass
        self.SMA = SMA
        self.J2 = J2
        self.GM = GM


#--------------OTHER STUFF --------------------------

# Magnetorquers


# Sensors (Noise constants?)

# class Moon():
#     """
#     A class to store Moon parameters
#     """
#     def __init__(self,
#                 R = 1738, #Equatorial Radius, km
#                 mass = 7.3483e22, #kg
#                 SMA = 38400, #semimajor axis, km
#                 GM = 4902.799): #gravitational param, km^3/s^2            
#         self.R = R
#         self.mass = mass
#         self.SMA = SMA
#         self.GM = GM
