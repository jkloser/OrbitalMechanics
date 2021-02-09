from math import *
import numpy as np
from class_orbit import orbit
import central_bodies

"""
### Enter initial radius and velocity
rx = float(input('Enter rx value: ').replace('\U00002013', '-'))
ry = float(input('Enter ry value: ').replace('\U00002013', '-'))
rz = float(input('Enter rz value: ').replace('\U00002013', '-'))
vx = float(input('Enter vx value: ').replace('\U00002013', '-'))
vy = float(input('Enter vy value: ').replace('\U00002013', '-'))
vz = float(input('Enter vz value: ').replace('\U00002013', '-'))
r0 = np.array([rx, ry, rz])
v0 = np.array([vx, vy, vz])
"""

r0 = np.array([0, 0, 1.0])
v0 = np.array([1.0, 0, 0])


### Choose a central body and units: sun_AU, sun_metric, earth_DU, earth_metric

satellite = orbit(central_bodies.earth_DU)

satellite.init_cond(r0, v0)
satellite.rv2elem()