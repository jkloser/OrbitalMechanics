"""
John Kloser
AENG 4150 Orbital Mechanics
class_orbit.py
Project 1
"""

import numpy as np
import central_bodies
#import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

### Add a plot method

class orbit:
    def __init__(self, r, v, cb = central_bodies.earth_DU, tol = 5):
        # Initialize orbit variables used to calculate orbital elements
        self.r0 = r
        self.v0 = v
        self.cb = cb
        self.tol = tol

    def rv2elem(self, print_val=True):
        # Calculate specific energy, specific angular momentum, eccentricity vector, line of nodes
        energy = np.linalg.norm(self.v0)**2/2 - self.cb['mu']/np.linalg.norm(self.r0)
        h = np.cross(self.r0, self.v0)
        e = np.cross(self.v0, h)/self.cb['mu'] - self.r0/np.linalg.norm(self.r0)
        n = np.cross([0,0,1], h)

        # Calculate the semi-major axis. Returns 0 for a parabolic orbit.
        if round(energy,self.tol) == 0:
            a = float('NaN')
        else:
            a = -self.cb['mu']/(2*energy)

        # Calculate eccentricity magnitude
        e_mag = np.linalg.norm(e)

        # Calculate inclination of the orbit
        i = np.arccos(h[2]/np.linalg.norm(h))

        # Calculate argument of periapse
        # (angle between the line of nodes and eccentricity vector)
        if round(i,self.tol) == 0 or round(e_mag,self.tol) == 0 or round(i,self.tol) == 180: arg_periapse = float('NaN')
        else:
            arg_periapse = np.arccos(np.dot(n, e) / (np.linalg.norm(n)*e_mag))
            if e[2] < 0: arg_periapse = (2*np.pi) - arg_periapse    # quadrant check


        # Calculate longitude of ascending node 
        # (angle between I and line of nodes)
        if round(i,self.tol) == 0 or round(i,self.tol) == 180: long_AN = float('NaN')
        else:
            long_AN = np.arccos(n[0]/np.linalg.norm(n))
            if n[1] < 0: long_AN = (2*np.pi) - long_AN   # quadrant check

        # Calculate the orbit period
        if energy >= 0: period = float('NaN')
        else:
            period = 2*np.pi/np.sqrt(self.cb['mu']) * a**(3/2)

        # Calculate true anamoly of initial condiitons 
        # (angle between periapse and current position)
        if round(e_mag,self.tol) == 0: true_anamoly = float('NaN')
        else:
            true_anamoly = np.arccos(np.dot(e, self.r0) / (e_mag*np.linalg.norm(self.r0)))
            if np.dot(self.r0, self.v0) < 0: true_anamoly = (2*np.pi) - true_anamoly  # quadrant check

        # Calculate true longitude of periapse
        # (angle between I and eccentricity/periapse)
        if round(e_mag,self.tol) == 0: true_long_periapse = float('NaN')
        else:
            true_long_periapse = np.arccos(e[0]/np.linalg.norm(e))
            if e[1] < 0: true_long_periapse = (2*np.pi) - true_long_periapse     # quadrant check

        # Calculate the argument of latitude at epoch
        # (angle between line of nodes and satellite position)
        if round(i,self.tol) == 0 or round(i,self.tol) == 180: u = float('NaN')
        else:
            u = np.arccos(np.dot(n, self.r0) / (np.linalg.norm(n)*np.linalg.norm(self.r0)))
            if self.r0[2] < 0: u = (2*np.pi) - u  # quadrant check

        # Calculate the true longitude at epoch
        # (angle between I and satellite position)
        true_long_epoch = np.arccos(self.r0[0] / np.linalg.norm(self.r0))
        if self.r0[1] < 0: true_long_epoch = (2*np.pi) - true_long_epoch  # quadrant check

        # Calculate the eccentric anamoly and time since perigee
        k = np.sqrt(self.cb['mu']/a**3)
        if e_mag < 1:
            E = 2 * np.arctan(np.sqrt((1-e_mag)/(1+e_mag))*np.tan(true_anamoly/2))
            #E = E%(2*np.pi)
            if np.dot(self.r0, self.v0) < 0: E = 2*np.pi + E
            Tp = (e_mag*np.sin(E) - E)/k

        elif e_mag > 1:
            E = np.arctanh(np.sqrt((e_mag-1)/(e_mag+1)) * np.tan(true_anamoly/2))
            #E = E%(2*np.pi)
            # Double Check this quadrant check
            if np.dot(self.r0, self.v0) < 0: E = 2*np.pi + E  # quadrant check
            Tp = (e_mag*np.sin(E) - E)/k
        else:
            Tp = float('NaN')

        r2d = 180/np.pi

        if print_val:
            print('a = ' + str(a))
            print('e = ' + str(e_mag))
            print('i = ' + str(i*r2d))
            print('\u03C9 = ' + str(arg_periapse*r2d))
            print('\u03A9 = ' + str(long_AN*r2d))
            print('period = ' + str(period))
            print('\u03BD = ' + str(true_anamoly*r2d))
            print('\u03A0 = ' + str(true_long_periapse*r2d))
            print('u = ' + str(u*r2d))
            print('l = ' + str(true_long_epoch*r2d))
            print('Tp = ' + str(Tp))

        return np.array([a, e_mag, i, arg_periapse, long_AN, period, true_anamoly, true_long_periapse, u, true_long_epoch, Tp])

    def rotation_matrix(self):
        elem = self.rv2elem(print_val=False)

        i = elem[3]
        arg_periapse = elem[4]
        long_AN = elem[5]

        ### Special Cases
        # Circular Equatorial
        if elem[1] == 0 and elem[2] == 0:
            arg_periapse = 0
        # Circular Inclined
        elif elem[2] == 0:
            long_AN = 0
            arg_periapse = elem[7]
        # Elliptical Equatorial
        elif elem[1] == 0:
            arg_periapse = 0



        R11 = cos(long_AN)*cos(arg_periapse) - sin(long_AN)*sin(arg_periapse)*cos(i)
        R12 = -cos(long_AN)*sin(arg_periapse) - sin(long_AN)*cos(arg_periapse)*cos(i)
        R13 = sin(long_AN)*sin(i)

        R21 = sin(long_AN)*cos(arg_periapse) + cos(long_AN)*sin(arg_periapse)*cos(i)
        R22 = -sin(long_AN)*sin(arg_periapse) + cos(long_AN)*cos(arg_periapse)*cos(i)
        R23 = -cos(long_AN)*sin(i)

        R31 = sin(arg_periapse)*sin(i)
        R32 = cos(arg_periapse)*sin(i)
        R33 = cos(i)

        return np.array([[R11, R12, R13],
                        [R21, R22, R23],
                        [R31, R32, R33]])

def vehicle2rv(r, v, phi, Az, delta, GMST, lambdaE = 0, d2r = False):

    if d2r == True:
        conv = np.pi/180
        phi *= conv
        Az *= conv
        delta *= conv
        GMST *= conv
        lambdaE *= conv

    # Transforms vehicle centered 
    lst = GMST+lambdaE
    rx = r*np.cos(delta)*np.cos(lst)
    ry = r*np.cos(delta)*np.sin(lst)
    rz = r*np.sin(delta)

    v_S = -v*np.cos(phi)*np.cos(Az)
    v_E = v*np.cos(phi)*np.sin(Az)
    v_R = v*np.sin(phi)

    vx = v_S*np.sin(delta)*np.cos(lst) - v_E*np.sin(lst) + v_R*np.cos(delta)*np.cos(lst)
    vy = v_S*np.sin(delta)*np.sin(lst) + v_E*np.cos(lst) + v_R*np.cos(delta)*np.sin(lst)
    vz = -v_S*np.cos(delta) + v_R*np.sin(delta)
    
    return np.array([rx, ry, rz]), np.array([vx, vy, vz])


