from math import *
import numpy as np
from class_orbit import orbit

r0 = np.array([1, 2, 3])
v0 = np.array([3, 2, 1])
mu = 1

satellite = orbit(mu)
print(r0)
print(v0)
satellite.init_cond(r0, v0)
