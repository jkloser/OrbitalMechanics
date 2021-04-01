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
rx = float(input('Enter rx value: ').replace('\U00002013', '-'))
ry = float(input('Enter ry value: ').replace('\U00002013', '-'))
rz = float(input('Enter rz value: ').replace('\U00002013', '-'))
vx = float(input('Enter vx value: ').replace('\U00002013', '-'))
vy = float(input('Enter vy value: ').replace('\U00002013', '-'))
vz = float(input('Enter vz value: ').replace('\U00002013', '-'))
r0 = np.array([rx, ry, rz])
v0 = np.array([vx, vy, vz])

cb_choose = input('Choose earth_DU, earth_metric, sun_AU, sun_metric: ')
cb_list = {'earth_DU':cb.earth_DU, 'earth_metric':cb.earth_metric, 'sun_AU':cb.sun_AU, 'sun_metric':cb.sun_metric}


# r0 = np.array([-0.2851, 0.9411, 0.0])
# v0 = np.array([-0.9735, -0.2936, 0.0])
# r0 = np.array([1,0,0])
# v0 = np.array([0,2,0])
# cb_choose = 'earth_DU'

# DU = 6378.145
# DUTU = 7.905368

satellite = orbit(r0, v0, cb_list[cb_choose])

satellite.rv2elem()

dt = float(input('Enter time step: ').replace('\U00002013', '-'))
# dt = 1.7202 # TU
[r1, v1] = satellite.univ_formulation(dt)