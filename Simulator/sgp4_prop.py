# -*- coding: utf-8 -*-
import numpy as np
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv
import pdb

# Example ISS (Zarya) TLE
line1 = ('1 25544U 98067A   19302.69799672  .00001237  00000-0  29468-4 0  9996')
line2 = ('2 25544  51.6449  54.3108 0006382 201.3316 303.7490 15.50231045196128')

ISS = twoline2rv(line1, line2, wgs84)
position, velocity = ISS.propagate(2019, 10, 30, 00, 00, 00)

print(position)
print(velocity)