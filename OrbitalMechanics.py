from math import *
import numpy as np
from class_orbit import orbit
import central_bodies

### Enter initial radius and velocity
r0 = np.array([2, 0, 0])
v0 = np.array([0, 1, 0])

### Choose a central body and units: sun_AU, sun_metric, earth_DU, earth_metric

satellite = orbit(central_bodies.earth_DU)

satellite.init_cond(r0, v0)

print('a = ' + str(satellite.calc_a()))
print('e = ' + str(satellite.calc_e()))
print('i = ' + str(satellite.calc_i()))
print('\u03C9 = ' + str(satellite.calc_arg_periapse()))
print('\u03A9 = ' + str(satellite.calc_long_AN()))
print('period = ' + str(satellite.calc_time_period()))
print('\u03BD = ' + str(satellite.calc_true_anamoly()))
print('\u03A0 = ' + str(satellite.calc_true_long_periapse()))
print('u = ' + str(satellite.calc_arg_lat_epoch()))
print('l = ' + str(satellite.calc_true_long_epoch()))
