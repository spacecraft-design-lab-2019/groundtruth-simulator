from astropy_sun_position import sun_position_ECI
from astropy_sun_position import approx_sun_position_ECI
import math
import numpy as np

# angle between two vectors

def test_sun_model():

	MJD = 54000
	a = approx_sun_position_ECI(MJD)
	b = sun_position_ECI(MJD)
	angle_diff = math.acos(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)))*180/math.pi

	assert angle_diff < 5
