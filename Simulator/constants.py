# -*- coding: utf-8 -*-
import numpy as np
import msise00
import julian
import conversions as conv
import math
import pyIGRF

#--------------------Spacecraft Structure----------------------------

class SpacecraftStructure():
    """
    A class to store spacecraft structural poperties
    """
    def __init__(self,
                 I = np.array([[17,0,0],[0,18,0],[0,0,22]]),
                 cD = 2.3): #drag coefficient
        self.I = I
        self.cD = cD
        self.faces = self.make_faces()
        # self.surfArea = sum(f.A for f in faces)

    def aerodrag(self, rho, vRel):
        """
        Return drag force (F) and moment (M) due to atmospheric drag
         - rho is the local atmospheric density.
         - vRel is velocity relative to the atmosphere. Must be in the body frame.
        """
        F = np.zeros(3)
        M = np.zeros(3)

        for face in self.faces:
            f, m = face.aerodrag(rho, vRel)
            F += f
            M += m
        return self.cD*F, self.cD*M

    def make_faces(self):
        L = 0.05 # 50 mm sidelength
        c_len = L/2 # distance from center to a face
        A_main = L**2
        A_long = 0.064 * 0.058 # area of the special face

        X = unit('x')
        Y = unit('y')
        Z = unit('z')
        return [
                Face(X,  A_main, c_len*X),
                Face(Y,  A_main, c_len*Y),
                Face(Z,  A_main, c_len*Z),
                Face(-X, A_main, -c_len*X),
                Face(-Y, A_long, -c_len*Y),
                Face(-Z, A_main, -c_len*Z),
                # special annular faces. NOTE: these only face +y. The -y is taken care of by an extra-big face at (5)
                Face(Y, 0.004*0.064, np.array([0.027, -0.025, 0])),
                Face(Y, 0.004*0.064, np.array([-0.027, -0.025, 0])),
                Face(Y, 0.007*0.050, np.array([0, -0.025, 0.0285])),
                Face(Y, 0.007*0.050, np.array([0, -0.025, -0.0285]))
                ]


class Face():
    def __init__(self, N, A, c):
        self.N = N/np.linalg.norm(N)
        self.A = A
        self.c = c # vector from COM to CP of face.

    def wetted_area(self, v):
        """
        Compute the effective area of a face relative to a vector v where
        v is defined to point outwards from the face.
        """
        a = (v/np.linalg.norm(v)) @ self.N
        return max(a, 0) * self.A

    def aerodrag(self, rho, vRel):
        """
        See SpacecraftStruct.aerodrag()
        """
        vmag = np.linalg.norm(vRel)
        A = self.wetted_area(vRel)

        drag_acc = -0.5*rho*A*vmag * vRel
        drag_M = np.cross(self.c, drag_acc);
        return drag_acc, drag_M


#-------------------------Environment---------------------------------

class Environment():
    """
    A class to store environment constants / lookup functions.
    """
    def __init__(self, mjd_start = 58777.740671):
        self.mjd_start = mjd_start
        self.earth = Earth()

    def density_lookup(self, r_ECI, GMST, mjd):
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
        r_ECEF = conv.ECI_to_ECEF(r_ECI, GMST)
        glat, glon, alt = conv.ECEF_to_LLA(r_ECEF, self.earth.radius)
        t = julian.from_jd(mjd, fmt='mjd')

        atmos = msise00.run(time=t, altkm=alt, glat=glat, glon=glon)
        rho = atmos.Total.values[0].item()
        return rho

    def magfield_lookup(self, r_ECI, GMST, mjd):
        """
        Function: magfield_lookup

        Gets magnetic field vector using IGRF-12 model.
        Must have https://pypi.org/project/pyIGRF/ library installed.

        Inputs:
            year
            altitude (km)
            glat: geodetic latitude
            glon: geodetic longitude
        """
        r_ECEF = conv.ECI_to_ECEF(r_ECI, GMST)
        glat, glon, alt = conv.ECEF_to_LLA(r_ECEF, self.earth.radius)
        year = julian.from_jd(mjd, fmt='mjd').year

        field = pyIGRF.igrf_value(glat, glon, alt, year)
        B_NED = np.array([field[3], field[4], field[5]])
        return B_NED


class Earth():
    """
    A class to store Earth parameters
    """
    def __init__(self,
                radius = 6378.1378366, # Earth Equatorial Radius
                w = np.array([0, 0, 2*np.pi/86164.1]), # rotation rate vector of earth (in J2000, rotation is aligned with z-axis)
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

def unit(ax):
    """
    make a unit vector from int or str input (0, 1, 2) / ('x', 'y', 'z')
    """
    if isinstance(ax, str):
        if ax == 'x':
            ax = 0
        elif ax == 'y':
            ax = 1
        elif ax == 'z':
            ax = 2
        else:
            raise(Exception("invalid unit axis input {}".format(ax)))

    N = np.zeros(3)
    N[ax] = 1
    return N