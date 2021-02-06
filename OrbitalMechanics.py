from math import *
import numpy as np
from class_orbit import orbit
import central_bodies

### Enter initial radius and velocity
r0 = np.array([1, 2, 3])
v0 = np.array([3, 2, 1])

### Choose a central body and units: sun_AU, sun_metric, earth_DU, earth_metric

satellite = orbit(central_bodies.earth_DU)

satellite.init_cond(r0, v0)
