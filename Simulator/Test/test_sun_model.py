import pytest
import os, sys, inspect
import math
import numpy as np
import julian

# add current folder to the path
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
groundtruth_dir = os.path.dirname(currentdir)
sys.path.insert(0, groundtruth_dir)

import sun_model
from conversions import norm2

# angle between two vectors

def test_sun_position_ECI():

    start_date = 58890 # feb. 11, 2020 in MJD
    for i in range(0, 365, 5): # every 5 days for one year from start date
        mjd = start_date + i
        dt = julian.from_jd(mjd, fmt='mjd')

        sun_approx = sun_model.approx_sun_position_ECI(dt)
        sun_lookup = sun_model.sun_position_ECI(dt)
        angle_diff = math.acos(np.dot(sun_approx, sun_lookup)/(norm2(sun_approx)*norm2(sun_lookup)))*180/math.pi

        print(angle_diff)
        assert angle_diff < 5  # never exceed 5 degree error