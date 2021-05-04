"""
John Kloser
AENG 4150 Orbital Mechanics
OrbitalMechanics.py
Project 1
"""

import numpy as np
from class_orbit import *
import central_bodies as cb

### Enter initial radius and velocity. Replace the dash with a negative sign

r1x = float(input('Enter r1x value: ').replace('\U00002013', '-'))
r1y = float(input('Enter r1y value: ').replace('\U00002013', '-'))
r1z = float(input('Enter r1z value: ').replace('\U00002013', '-'))
r2x = float(input('Enter r2x value: ').replace('\U00002013', '-'))
r2y = float(input('Enter r2y value: ').replace('\U00002013', '-'))
r2z = float(input('Enter r2z value: ').replace('\U00002013', '-'))
"""
r1 = np.array([0.5, 0.6, 0.7])
r2 = np.array([0.0, 1.0, 0.0])
dt = 0.9667663
direction = -1
"""

r1 = np.array([r1x, r1y, r1z])
r2 = np.array([r2x, r2y, r2z])
dt = float(input('Enter time step: '))
direction = float(input('Enter 1 for short, -1 for long direction: ').replace('\U00002013', '-'))

cb_choose = input('Choose earth_DU, earth_metric, sun_AU, sun_metric: ')
cb_list = {'earth_DU':cb.earth_DU, 'earth_metric':cb.earth_metric, 'sun_AU':cb.sun_AU, 'sun_metric':cb.sun_metric}

[v1, v2] = gauss_problem(r1, r2, dt, direction, cb_list[cb_choose])
print('v1 = ' + str(v1))
print('v2 = ' + str(v2))