import numpy as np
import central_bodies as cb
from class_orbit import *

def rotation_matrix(i, long_AN, arg_periapse):
	i *= np.pi/180
	long_AN *- np.pi/180
	arg_periapse *= np.pi/180

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
JD_earth = 367.*year - np.fix(7.*(year+np.floor((month+9.)/12.))/4.) + np.fix(275.*month/9.) + day + 1721013.5 + ((sec/60. + minute)/60. + hour)/24.

# Arrival Daet at Venus (injection burn)
year = 2006.0
month = 04.0
day = 11.0
hour = 07.0
minute = 10.0
sec = 29.0
JD_venus = 367.*year - np.fix(7.*(year+np.floor((month+9.)/12.))/4.) + np.fix(275.*month/9.) + day + 1721013.5 + ((sec/60. + minute)/60. + hour)/24.

TOF = JD_venus - JD_earth # give time of flight in JD, need to conver to TU
#TOF /= 58.132821

T_earth = (JD_earth - 2451545.0)/36525.0
T_venus = (JD_venus - 2451545.0)/36525.0

a_earth = 1.000001018
e_earth = 0.01670862 - 0.000042037*T_earth
i_earth = 0.0000000 + 0.0130546*T_earth
Omega_earth = 0
omega_tilde_earth = 102.937348 + 0.3225557*T_earth
lambda_M_earth = 100.466449 + 35999.3728519*T_earth
omega_earth = omega_tilde_earth - Omega_earth

a_venus = 0.723329820
e_venus = 0.00677188- 0.000047788*T_venus
i_venus = 3.394662 - 0.0008568*T_venus
Omega_venus = 76.679920 - 0.278008*T_venus
omega_tilde_venus = 131.563707 + 0.0048646*T_venus
lambda_M_venus = 181.979801 + 58517.815676*T_venus
omega_venus = omega_tilde_venus - Omega_venus

lambda_M_earth %= 360
lambda_M_venus %= 360

M_earth = lambda_M_earth - omega_tilde_earth
M_venus = lambda_M_venus - omega_tilde_venus

M_earth *= np.pi/180
M_venus *= np.pi/180

M_earth %= 2*np.pi
M_venus %= 2*np.pi

E_earth = M_earth
iter = 1
print('E_earth')
while ((np.abs(E_earth - e_earth*np.sin(E_earth) - M_earth)) > 1*10**-9):
	E_earth = E_earth - (E_earth - e_earth*np.sin(E_earth) - M_earth) / (1-e_earth*np.cos(E_earth))
	iter += 1;
	print(iter)
	print(E_earth)

print('')
print('E_venus')
E_venus=M_venus
iter = 1
while ((np.abs(E_venus - e_venus*np.sin(E_venus) - M_venus)) > 1*10**-9):
	E_venus = E_venus - (E_venus - e_venus*np.sin(E_venus) - M_venus) / (1-e_venus*np.cos(E_venus))
	iter += 1
	print(iter)
	print(E_venus)

nu_earth = np.arccos(e_earth - np.cos(E_earth) ) / (e_earth*np.cos(E_earth)-1)
if ((E_earth > np.pi) & (nu_earth < np.pi)):
	nu_earth = 2*np.pi - nu_earth

nu_venus = np.arccos(e_venus - np.cos(E_venus)) / (e_venus*np.cos(E_venus)-1)
if ((E_venus > np.pi) & (nu_venus < np.pi)):
	nu_venus = 2*np.pi - nu_venus

p_earth = a_earth*(1 - e_earth**2)
r_earth_pqw = np.array([p_earth*np.cos(nu_earth)/(1+e_earth*np.cos(nu_earth)), p_earth*np.sin(nu_earth)/(1+e_earth*np.cos(nu_earth)), 0.0])
v_earth_pqw = np.sqrt(mu/p_earth) * np.array([-np.sin(nu_earth), e_earth+np.cos(nu_earth), 0.0])

r_earth_ijk = rotation_matrix(i_earth, Omega_earth, omega_earth) @ np.transpose(r_earth_pqw)
v_earth_ijk = rotation_matrix(i_earth, Omega_earth, omega_earth) @ np.transpose(v_earth_pqw)

p_venus = a_venus*(1 - e_venus**2)
r_venus_pqw = np.array([p_venus*np.cos(nu_venus)/(1+e_venus*np.cos(nu_venus)), p_venus*np.sin(nu_venus)/(1+e_venus*np.cos(nu_venus)), 0.0])
v_venus_pqw = np.sqrt(mu/p_venus) * np.array([-np.sin(nu_venus), e_venus+np.cos(nu_venus), 0.0])


r_venus_ijk = rotation_matrix(i_venus, Omega_venus, omega_venus) @ np.transpose(r_venus_pqw)
v_venus_ijk	= rotation_matrix(i_venus, Omega_venus, omega_venus) @ np.transpose(v_venus_pqw)

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
# check pqw is same as ijk
earth = orbit(np.transpose(r_earth_ijk), np.transpose(v_earth_ijk), cb.sun_AU)
venus = orbit(np.transpose(r_venus_ijk), np.transpose(v_venus_ijk), cb.sun_AU)

print('')
print('Earth:')
earth.rv2elem()

print('')
print('Venus:')
venus.rv2elem()