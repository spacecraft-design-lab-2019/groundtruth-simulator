from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body, get_moon, get_sun
import numpy as np

def sun_position_ECI(MJD):
	# this function gets the ECI sun position
	# INPUTS:
	# t in MJD
	# OUTPUTS:
	# r_ECI (X, Y, X) in (km, km, km)
	AU_km = 149597870.7

	t = Time(MJD, scale = 'utc',format = 'mjd')

	try:
		sun = get_sun(t)
	except NameError:
		print('Astropy not installed ya dum dum')

	r_GCRS = np.array([sun.ra.value, sun.dec.value, sun.distance.value * AU_km])

	r_ECI = GCRS_to_ECI(r_GCRS)

	return (r_ECI)

def GCRS_to_ECI(r_GCRS):
	# this function converts from GCRS to ECEF:
	# INPUTS:
	# r_GCRS ( ra, dec, distance) in (deg, deg, km)
	# OUTPUTS:
	# r_ECI (X, Y, Z) in (km, km, km)

	ra = r_GCRS[0] * np.pi / 180
	dec = r_GCRS[1] * np.pi / 180
	dist = r_GCRS[2]

	r_ECI = np.array([dist * np.cos(ra),
		 dist * np.sin(ra) * np.cos(dec),
		 dist * np.sin(ra) * np.sin(dec)])

	return (r_ECI)

t = Time(51622, scale = 'utc',format = 'mjd')

AU_km = 149597870.7

sun = get_sun(t)

r_GCRS = np.array([sun.ra.value, sun.dec.value, sun.distance.value * AU_km])

r_ECI = GCRS_to_ECI(r_GCRS)
