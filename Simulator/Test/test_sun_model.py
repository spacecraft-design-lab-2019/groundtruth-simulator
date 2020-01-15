import math
import numpy as np
import sun_model
import julian

# angle between two vectors

def test_sun_position_ECI():

	mjd = 54000
	dt = julian.from_jd(mjd, fmt='mjd')

	a = sun_model.approx_sun_position_ECI(mjd)
	b = sun_model.sun_position_ECI(dt)
	angle_diff = math.acos(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)))*180/math.pi

	assert angle_diff < 5