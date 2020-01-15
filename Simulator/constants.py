# -*- coding: utf-8 -*-
import numpy as np
# import msise00
import julian
import conversions as conv
import pyIGRF
import sun_model
import math
import sys
sys.path.append('/home/eleboeuf/Documents/GNC')
sys.path.append('./GNC')
import magnetic_field_cpp as mfcpp
#--------------------Spacecraft Structure----------------------------

class SpacecraftStructure():
    """
    A class to store spacecraft structural poperties
    """
    def __init__(self, I, cD=2.3, mass=1.0):
        self.I = I
        self.cD = cD
        self.mass = mass # kg
        self.faces = self.make_faces()

    def aerodrag(self, rho, vRel):
        """
        Return drag force (F) and moment (M) due to atmospheric drag
         - rho is the local atmospheric density.
         - vRel is velocity relative to the atmosphere. Must be in the body frame.
        """
        F = np.zeros(3)
        M = np.zeros(3)
        # vmag = np.linalg.norm(vRel)
        vmag = (vRel.T @ vRel)**(0.5)
        for face in self.faces:
            f, m = face.aerodrag(rho, vRel, vmag)
            F += f
            M += m
        return self.cD/self.mass*F, self.cD/self.mass*M

    def make_faces(self):
        L = 0.05 # 50 mm sidelength
        c_len = L/2 # distance from center to a face
        A_main = L**2
        A_long = 0.064 * 0.058 # area of the special face

        X = conv.unit('x')
        Y = conv.unit('y')
        Z = conv.unit('z')
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

    def wetted_area(self, v_unit):
        """
        Compute the effective area of a face relative to a vector v where
        v is the unit vector defined to point outwards from the face.
        """
        # a = (v/np.linalg.norm(v)) @ self.N
        a = v_unit @ self.N
        return max(a, 0) * self.A

    def aerodrag(self, rho, vRel, vmag):
        """
        See SpacecraftStruct.aerodrag()
        """
        # vmag = np.linalg.norm(vRel)
        A = self.wetted_area(vRel/vmag)

        drag_acc = -0.5*rho*A*vmag * vRel
        # drag_M = np.cross(self.c, drag_acc);

        drag_M = conv.cross3(self.c, drag_acc)
        return drag_acc, drag_M


#-------------------------Environment---------------------------------

class Environment():
    """
    A class to store environment constants / lookup functions.
    """
    def __init__(self, datetime):
        self.datetime = datetime
        self.earth = Earth()

    def update(self, datetime_new):
        self.datetime = datetime_new

    def density_lookup(self, r_ECI, model="exponential_2"):
        """
        Function: density_lookup

        Gets atmospheric density using various models.
        For MSISE-00 atmospheric model, must have https://pypi.org/project/msise00/ library installed.

        Inputs:
            r_ECI:  position vector in ECI frame
            model:  parameter to select type of model (default="exponential")

        Outputs:
            rho: atmospheric density (kg/m^3)
        """

        if model == "exponential":
            rho_0 = 1.225e9
            h = np.linalg.norm(r_ECI)
            H = 10.0
            rho = rho_0 * np.exp(-(h-self.earth.radius)/H)

        elif model == "exponential_2":
            #drag equation fit coefficients
            a = 4.436e-09
            b = -0.01895
            c = 4.895e-12
            d = -0.008471
            # R = np.linalg.norm(r_ECI) - self.earth.radius
            R = conv.norm2(r_ECI) - self.earth.radius
            rho = a*np.exp(b*R) + c*np.exp(d*R)

        elif model == "msise00":
            mjd = julian.to_jd(self.datetime, fmt='mjd')
            GMST = conv.mjd_2_GMST(mjd)
            r_ECEF = conv.ECI_to_ECEF(r_ECI, GMST)
            glat, glon, alt = conv.ECEF_to_LLA(r_ECEF, self.earth.radius)

            atmos = msise00.run(time=self.datetime, altkm=alt, glat=glat, glon=glon)
            rho = atmos.Total.values[0].item()

        return rho

    def magfield_lookup(self, r_ECI, order):
        """
        Function: magfield_lookup

        Gets magnetic field vector using IGRF-12 model.
        Must have https://pypi.org/project/pyIGRF/ library installed.

        Inputs:
            year
            altitude (km)
            glat: geodetic latitude
            glon: geodetic longitude
        Ouputs:
            B_NED: magnetic field vector in North/East/Down
        """
        mjd = julian.to_jd(self.datetime, fmt='mjd')
        GMST = conv.mjd_2_GMST(mjd)
        r_ECEF = conv.ECI_to_ECEF(r_ECI, GMST)
        glat, glon, alt = conv.ECEF_to_LLA(r_ECEF, self.earth.radius)
        lat = np.rad2deg(glat)
        lon = np.rad2deg(glon) % 360.0                      # pyIGRF requires East Longitude (0-360 deg)
        year = self.datetime.year

        # field = pyIGRF.igrf_value(lat, lon, alt, year)      # pyIGRF uses degrees!!!!

        B_NED = mfcpp.get_magnetic_field(lat, lon, alt, year, order)

        # B_NED = np.array([field[3], field[4], field[5]])
        B_ECI = conv.NED_to_ECI(B_NED, glat, glon, GMST)
        return B_ECI

    def sunVector(self, r_ECI):
        """
        Inputs:
            r_ECI: ECI position of satellite
        Output:
            position vector (3-vector) from satellite to sun
        """
        jdate = julian.to_jd(self.datetime, fmt='jd')
        earth2sun = sun_model.sun1(jdate)
        rsun = (earth2sun - r_ECI)/np.linalg.norm(earth2sun - r_ECI)
        return rsun


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
