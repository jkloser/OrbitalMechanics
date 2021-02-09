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

r0 = np.array([-6045, -3490, 2500])
v0 = np.array([-3.457, 6.618, 2.533])

DU = 6378.145
DUTU = 7.905368


### Choose a central body and units: sun_AU, sun_metric, earth_DU, earth_metric

satellite = orbit(r0/DU, v0/DUTU, central_bodies.earth_DU)

satellite.rv2elem()