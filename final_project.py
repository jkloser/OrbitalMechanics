import numpy as np
import central_bodies as cb
from class_orbit import *

def rotation_matrix(i, long_AN, arg_periapse):
	R11 = np.cos(long_AN)*np.cos(arg_periapse) - np.sin(long_AN)*np.sin(arg_periapse)*np.cos(i)
	R12 = -np.cos(long_AN)*np.sin(arg_periapse) - np.sin(long_AN)*np.cos(arg_periapse)*np.cos(i)
	R13 = np.sin(long_AN)*np.sin(i)

	R21 = np.sin(long_AN)*np.cos(arg_periapse) + np.cos(long_AN)*np.sin(arg_periapse)*np.cos(i)
	R22 = -np.sin(long_AN)*np.sin(arg_periapse) + np.cos(long_AN)*np.cos(arg_periapse)*np.cos(i)
	R23 = -np.cos(long_AN)*np.sin(i)

	R31 = np.sin(arg_periapse)*np.sin(i)
	R32 = np.cos(arg_periapse)*np.sin(i)
	R33 = np.cos(i)

	return np.array([[R11, R12, R13],
		[R21, R22, R23],
        [R31, R32, R33]])

# Departure Date from Earth
mu = 1.0
year = 2005.0
month = 11.0 
day = 09.0
hour = 05.0
minute = 09.0
sec = 00.0
JD_dep = 367.*year - np.fix(7.*(year+np.floor((month+9.)/12.))/4.) + np.fix(275.*month/9.) + day + 1721013.5 + ((sec/60. + minute)/60. + hour)/24.

# Arrival Daet at Venus (injection burn)
year = 2006.0
month = 04.0
day = 11.0
hour = 07.0
minute = 10.0
sec = 29.0
JD_arr = 367.*year - np.fix(7.*(year+np.floor((month+9.)/12.))/4.) + np.fix(275.*month/9.) + day + 1721013.5 + ((sec/60. + minute)/60. + hour)/24.

TOF = JD_arr - JD_dep

TOF = JD_arr - JD_dep

T_dep = (JD_dep - 2451545.0)/36525.0
T_arr = (JD_arr - 2451545.0)/36525.0

a_dep = 1.000001018
e_dep = 0.01670862 - 0.000042037*T_dep - 0.0000001236*T_dep**2 + 0.00000000004*T_dep**3
i_dep = 0.0000000 + 0.0130546*T_dep - 0.00000931*T_dep**2 - 0.000000034*T_dep**3
Omega_dep = 0
omega_dep = 102.937348 + 0.322555*T_dep + 0.00015026*T_dep**2 + 0.000000478*T_dep**3
lambda_M_dep = 100.466449 + 35999.3728519*T_dep - 0.00000568*T_dep**2

a_arr = 0.723329820
e_arr = 0.00677188- 0.000047788*T_arr + 0.0000000975*T_arr**2 + 0.00000000044*T_arr**3
i_arr = 3.394662 - 0.0008568*T_arr - 0.00003244*T_arr**2 + 0.000000010*T_arr**3
Omega_arr = 76.679920 - 0.278008*T_arr - 0.00014256*T_arr**2 - 0.000000198*T_arr**3
omega_arr = 131.563707 + 0.0048646*T_arr - 0.00138232*T_arr**2 - 0.000005332*T_arr**3
lambda_M_arr = 181.979801 + 58517.815676*T_arr + 0.00000165*T_arr**2 - 0.000000002*T_arr**3

lambda_M_dep %= 360
lambda_M_arr %= 360

lambda_M_dep *= np.pi/180
lambda_M_arr *= np.pi/180

M_earth = lambda_M_dep - omega_dep
M_venus = lambda_M_arr - omega_arr

E_earth = M_earth
iter = 1
print('E_earth')
while ((np.abs(E_earth - e_dep*np.sin(E_earth) - M_earth)) > 1*10**-9):
	E_earth = E_earth - (E_earth - e_dep*np.sin(E_earth) - M_earth) / (1-e_dep*np.cos(E_earth))
	iter += 1;
	print(iter)
	print(E_earth)

print('')
print('E_venus')
E_venus=M_venus
iter = 1
while ((np.abs(E_venus - e_arr*np.sin(E_venus) - M_venus)) > 1*10**-9):
	E_venus = E_venus - (E_venus - e_arr*np.sin(E_venus) - M_venus) / (1-e_arr*np.cos(E_venus))
	iter += 1
	print(iter)
	print(E_venus)

nu_earth = np.arccos(np.cos(E_earth) - e_dep) / (1 - e_dep*np.cos(E_earth))
if ((E_earth > np.pi) & (nu_earth < np.pi)):
	nu_earth += np.pi

nu_venus = np.arccos(np.cos(E_venus) - e_arr) / (1 - e_arr*np.cos(E_venus))
if ((E_venus > np.pi) & (nu_venus < np.pi)):
	nu_venus += np.pi

p_earth = a_dep*(1 - e_dep**2)
r_earth_pqw = np.array([p_earth*np.cos(nu_earth)/(1+e_dep*np.cos(nu_earth)), p_earth*np.sin(nu_earth)/(1+e_dep*np.cos(nu_earth)), 0.0])
v_earth_pqw = np.sqrt(mu/p_earth) * np.array([-np.sin(nu_earth), e_dep+np.cos(nu_earth), 0.0])

r_earth_ijk = rotation_matrix(i_dep, Omega_dep, omega_dep) @ np.transpose(r_earth_pqw)
v_earth_ijk = rotation_matrix(i_dep, Omega_dep, omega_dep) @ np.transpose(v_earth_pqw)

p_venus = a_arr*(1 - e_arr**2)
r_venus_pqw = np.array([p_venus*np.cos(nu_venus)/(1+e_arr*np.cos(nu_venus)), p_venus*np.sin(nu_venus)/(1+e_arr*np.cos(nu_venus)), 0.0])
v_venus_pqw = np.sqrt(mu/p_venus) * np.array([-np.sin(nu_venus), e_arr+np.cos(nu_venus), 0.0])

r_venus_ijk = rotation_matrix(i_arr, Omega_arr, omega_arr) @ np.transpose(r_venus_pqw)
v_venus_ijk	= rotation_matrix(i_arr, Omega_arr, omega_arr) @ np.transpose(v_venus_pqw)

print('')
print('Earth')
print(r_earth_ijk)
print(v_earth_ijk)

print('')
print('Venus')
print(r_venus_ijk)
print(v_venus_ijk)
"""
r_dep = np.array([0.969829190018906, 0.278019514681994, 0.000212511971878])
r_arr = np.array([-0.564825421309492, -0.439595103188768, 0.104932439376776])

v_dep = np.array([-0.286578530958222, 0.948812682284867, 0.000725251442459])
v_arr = np.array([-0.540379616661920, -1.038183441347880, 0.052202483476656])
"""

earth = orbit(np.transpose(r_earth_ijk), np.transpose(v_earth_ijk), cb.sun_AU)
venus = orbit(np.transpose(r_venus_ijk), np.transpose(v_venus_ijk), cb.sun_AU)

print('')
print('Earth:')
earth.rv2elem()

print('')
print('Venus:')
venus.rv2elem()