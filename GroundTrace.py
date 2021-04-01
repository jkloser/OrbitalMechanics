import numpy as np
from class_orbit import *
import central_bodies as cb

r = 4263.195*5280
v = 25000
phi = 3
lambdaE = 279.45
delta = 28.5
Az = 75
gmst = 286.489961


r, v = vehicle2rv(r, v, phi, Az, delta, gmst, lambdaE, True)

print(r)
print(v)

# satellite = orbit(r, v, cb.earth_imperial)
# satellite.rv2elem()