# -*- coding: utf-8 -*-
import numpy as np
# import msise00
import julian
import conversions as conv
import pyIGRF
import sun_model
import math
import sys, os
import sun_sensor_math
from collections import namedtuple

dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, dir+'/GNC/')
import magnetic_field_cpp as mfcpp
#--------------------Spacecraft Structure----------------------------

class SpacecraftStructure():
    """
    A class to store spacecraft structural poperties
    """
    def __init__(self, I, mass, thermal_properties, cD=2.3):
        self.I = I
        self.cD = cD
        self.mass = mass
        self.thermal_properties = thermal_properties
        self.faces = self.make_faces()

    def aerodrag(self, rho, vRel):
        """
        Return acceleration due to drag (F/m) and moment-per-mass (M/m) due to atmospheric drag
         - rho is the local atmospheric density.
         - vRel is velocity relative to the atmosphere *in the body frame*.
        """
        F = np.zeros(3)
        M = np.zeros(3)
        vmag = conv.norm2(vRel)
        for face in self.faces:
            # If the velocity relative to the face normal is negative,
            # it is not acted on by drag and we can skip it.
            a = vRel @ face.N
            if a > 0:
                # Drag = -1/2 * Cd * ρ*A_eff*v² * v̂
                # constants that are the same across the faces are multiplied
                # at the end so we have: F = A_eff*v_vec
                f = face.A*a * vRel
                F += f
                M += conv.cross3(face.c, f)

        C = -0.5 * rho * vmag * self.cD / self.mass
        return C*F, C*M

    def make_faces(self):
        L = 0.05 # 50 mm sidelength
        c_len = L/2 # distance from center to a face
        A_main = L**2
        A_long = 0.064 * 0.058 # area of the special face

        X = conv.unit('x')
        Y = conv.unit('y')
        Z = conv.unit('z')

        # Named tuple to replace the Face class.
        Face = namedtuple('Face', 'N, A, c')

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

    def calc_Tdot(self, T, eclipse):
        alpha = self.thermal_properties["absorptivity"]
        eps = self.thermal_properties["emmisivity"]
        A = self.faces[0].A
        Js = 1368 # W / m^2
        Ja = Js * self.thermal_properties["albedo"]
        Je = 231 # W / m^2
        sigma = 5.67e-8; # Stefan Bolzmann constant

        if eclipse:
            Q_in = alpha * A * Je
        else:
            Q_in = alpha * A * (Js + Ja + Je)

        Tdot = (Q_in - eps*sigma*(T**4)*6*A) / (self.mass * self.thermal_properties["heat_capacitance"])
        return Tdot



#-------------------------Environment---------------------------------

class Environment():
    """
    A class to store environment constants / lookup functions.
    """
    def __init__(self, mjd_start):
        self.mjd = mjd_start
        self.mjd_start = mjd_start
        self.earth = Earth()

    def update(self, t):
        self.mjd = self.mjd_start + t/(60*60*24)

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
            GMST = conv.mjd_2_GMST(self.mjd)
            r_ECEF = conv.ECI_to_ECEF(r_ECI, GMST)
            glat, glon, alt = conv.ECEF_to_LLA(r_ECEF, self.earth.radius)

            datetime = julian.from_jd(self.mjd, fmt='mjd')
            atmos = msise00.run(time=datetime, altkm=alt, glat=glat, glon=glon)
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
            B_NED: magnetic field vector in North/East/Down in units of Teslas
        """
        GMST = conv.mjd_2_GMST(self.mjd)
        r_ECEF = conv.ECI_to_ECEF(r_ECI, GMST)
        glat, glon, alt = conv.ECEF_to_LLA(r_ECEF, self.earth.radius)
        lat = np.rad2deg(glat)
        lon = np.rad2deg(glon) % 360.0                      # pyIGRF requires East Longitude (0-360 deg)
        year = julian.from_jd(self.mjd, fmt='mjd').year

        # field = pyIGRF.igrf_value(lat, lon, alt, year)      # pyIGRF uses degrees!!!!

        B_NED = mfcpp.get_magnetic_field(lat, lon, alt, year, order)/1e9

        # B_NED = np.array([field[3], field[4], field[5]])
        B_ECI = conv.NED_to_ECI(B_NED, glat, glon, GMST)
        return B_ECI

    def sunVector(self, r_ECI):
        """
        Inputs:
            r_ECI: ECI position of satellite
        Output:
            position unit vector (3-vector) from satellite to sun
        """
        earth2sun = sun_model.approx_sun_position_ECI(self.mjd)
        rsun = (earth2sun - r_ECI)/np.linalg.norm(earth2sun - r_ECI)
        return rsun

    def isEclipse(self, r_ECI):
        """
        Determines if the satellite is in eclipse

        Inputs:
            r_ECI: ECI position of the satellite
        Output:
            boolean: (true if in eclipse)
        """
        rsun = self.sunVector(r_ECI)
        return sun_sensor_math.isEclipse(r_ECI, rsun, self.earth.radius)


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
