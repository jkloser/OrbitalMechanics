import numpy as np
import central_bodies as cb 
from class_orbit import *

r0 = np.array([-0.2851, 0.9411, 0.0])
v0 = np.array([-0.9735, -0.2936, 0.0])
cb_choose = 'earth_DU'

satellite = orbit(r0, v0, cb_list[cb_choose])

satellite.rv2elem()


dt = 1.7202 # TU
[r1, v1] = satellite.univ_formulation(dt)

# r=
# v=
# phi=
# Az=
# delta=
# GMST=
# [r_vec, v_vec] = vehicle2rv(r, v, phi, Az, delta, GMST, lambdaE=0, d2r=True)
# print('r = ' + str(r[0]) + 'i + ' str(r[1]) + 'j + ' + str(r[2]) + 'k')
# print('v = ' + str(v[0]) + 'i + ' str(v[1]) + 'j + ' + str(v[2]) + 'k')
