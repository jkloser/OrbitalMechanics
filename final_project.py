import numpy as np
import central_bodies as cb
from class_orbit import *

def rotation_matrix(i, long_AN, arg_periapse):
	i *= np.pi/180.0
	long_AN *- np.pi/180.0
	arg_periapse *= np.pi/180.0

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
year = 2018
month = 8
day = 12
hour = 7
minute = 31
sec = 00.0
JD_earth = 367.*year - np.fix(7.*(year+np.floor((month+9.)/12.))/4.) + np.fix(275.*month/9.) + day + 1721013.5 + ((sec/60. + minute)/60. + hour)/24.

# Arrival Daet at Venus (injection burn)
year = 2018
month = 10
day = 3
hour = 8
minute = 44
sec = 0.0
JD_venus = 367.*year - np.fix(7.*(year+np.floor((month+9.)/12.))/4.) + np.fix(275.*month/9.) + day + 1721013.5 + ((sec/60. + minute)/60. + hour)/24.

TOF = JD_venus - JD_earth # give time of flight in JD, need to conver to TU
dt = TOF / 58.132821

# Centuries since 2000
T_earth = (JD_earth - 2451545.0)/36525.0
T_venus = (JD_venus - 2451545.0)/36525.0

# Earth Ephemerides Ref: Meeus (1991, 202-204)
a_earth = 1.000001018
e_earth = 0.01670862 - 0.000042037*T_earth
i_earth = 0.0000000 + 0.0130546*T_earth
Omega_earth = 0
omega_tilde_earth = 102.937348 + 0.3225557*T_earth
lambda_M_earth = 100.466449 + 35999.3728519*T_earth
omega_earth = omega_tilde_earth - Omega_earth

# Venus Ephemerides Ref: Meeus (1991, 202-204)
a_venus = 0.723329820
e_venus = 0.00677188- 0.000047766*T_venus
i_venus = 3.394662 - 0.0008568*T_venus
Omega_venus = 76.679920 - 0.278008*T_venus
omega_tilde_venus = 131.563707 + 0.0048646*T_venus
lambda_M_venus = 181.979801 + 58517.815676*T_venus
omega_venus = omega_tilde_venus - Omega_venus

# Bring within 2pi range
lambda_M_earth %= 360
lambda_M_venus %= 360

# Find mean anomaly
M_earth = lambda_M_earth - omega_tilde_earth
M_venus = lambda_M_venus - omega_tilde_venus

M_earth *= np.pi/180
M_venus *= np.pi/180

M_earth %= 2*np.pi
M_venus %= 2*np.pi

# Iterate using Newton's Method to find eccentric anomaly
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

# Find true anomaly, checking to be in the same half plane as eccentric anomaly
nu_earth = np.arccos((np.cos(E_earth) - e_earth) / (1 - e_earth*np.cos(E_earth)))
#nu_earth = 2*np.arctan(np.tan(E_earth/2) / np.sqrt((1-e_earth)/(1+e_earth)))
if ((E_earth > np.pi) & (nu_earth < np.pi)):
	nu_earth = 2*np.pi - nu_earth

nu_venus = np.arccos((np.cos(E_venus) - e_venus) / (1 - e_venus*np.cos(E_earth)))
#nu_venus = 2*np.arctan(np.tan(E_venus/2) / np.sqrt((1-e_venus)/(1+e_venus)))
if ((E_venus > np.pi) & (nu_venus < np.pi)):
	nu_venus = 2*np.pi - nu_venus

# Find the earth position and velocity in perifocal coordinates
p_earth = a_earth*(1 - e_earth**2)
r_earth_pqw = np.array([p_earth*np.cos(nu_earth)/(1+e_earth*np.cos(nu_earth)), p_earth*np.sin(nu_earth)/(1+e_earth*np.cos(nu_earth)), 0.0])
v_earth_pqw = np.sqrt(mu/p_earth) * np.array([-np.sin(nu_earth), e_earth+np.cos(nu_earth), 0.0])

# Convert to heliocentric coordinates
r_earth_ijk = rotation_matrix(i_earth, Omega_earth, omega_earth) @ np.transpose(r_earth_pqw)
v_earth_ijk = rotation_matrix(i_earth, Omega_earth, omega_earth) @ np.transpose(v_earth_pqw)

# Venus perifocal vectors
p_venus = a_venus*(1 - e_venus**2)
r_venus_pqw = np.array([p_venus*np.cos(nu_venus)/(1+e_venus*np.cos(nu_venus)), p_venus*np.sin(nu_venus)/(1+e_venus*np.cos(nu_venus)), 0.0])
v_venus_pqw = np.sqrt(mu/p_venus) * np.array([-np.sin(nu_venus), e_venus+np.cos(nu_venus), 0.0])

rm = rotation_matrix(i_venus, Omega_venus, omega_venus)

r_venus_ijk = rm @ np.transpose(r_venus_pqw)
v_venus_ijk	= rotation_matrix(i_venus, Omega_venus, omega_venus) @ np.transpose(v_venus_pqw)

print('')
print('Earth')
print(r_earth_ijk)
print(v_earth_ijk)

print('')
print('Venus')
print(r_venus_ijk)
print(v_venus_ijk)

# Check the orbital parameters and compare to ephemerides
earth = orbit(np.transpose(r_earth_ijk), np.transpose(v_earth_ijk), cb.sun_AU)
venus = orbit(np.transpose(r_venus_ijk), np.transpose(v_venus_ijk), cb.sun_AU)

print('')
print('Earth:')
earth.rv2elem()

print('')
print('Venus:')
venus.rv2elem()

direction = -1
[transfer_orbit, v1, v2] = gauss_problem(r_venus_pqw, r_earth_pqw, dt, direction, cb.sun_AU)

transfer_orbit.rv2elem()
print('v1 = ' + str(v1))
print('v2 = ' + str(v2))