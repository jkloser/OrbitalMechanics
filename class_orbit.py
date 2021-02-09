import numpy as np
import central_bodies
import matplotlib as plt
### Add a plot method
### Check fringe cases: undefined angles, parabola, hyperbola, circular orbits
### tolerances for special cases (abs(i) < 0.0001)

class orbit:
    def __init__(self, cb = central_bodies.earth_DU):
        # Initialize orbit variables used to calculate orbital elements
        self.r0 = np.empty(3)
        self.v0 = np.empty(3)
        self.energy = np.empty(1)
        self.h = np.empty(3)
        self.e = np.empty(3)
        self.cb = cb

    def init_cond(self, r, v):
        # Set initial radius and velocity. Input as 1x3 numpy arrays.
        self.r0 = r
        self.v0 = v

    def rv2elem(self, print_val=True):
        # Calculate specific energy, specific angular momentum, eccentricity vector, line of nodes
        energy = np.linalg.norm(self.v0)**2/2 - self.cb['mu']/np.linalg.norm(self.r0)
        h = np.cross(self.r0, self.v0)
        e = np.cross(self.v0, h)/self.cb['mu'] - self.r0/np.linalg.norm(self.r0)
        n = np.cross([0,0,1], h)

        # Calculate the semi-major axis. Returns 0 for a parabolic orbit.
        if energy == 0:
            a = 0
        else:
            a = -self.cb['mu']/(2*energy)

        # Calculate eccentricity magnitude
        e_mag = np.linalg.norm(e)

        # Calculate inclination of the orbit
        i = np.arccos(h[2]/np.linalg.norm(h)) * 180 / np.pi

        # Calculate argument of periapse
        # (angle between the line of nodes and eccentricity vector)
        if i == 0 or e_mag == 0: arg_periapse = None
        else:
            arg_periapse = np.arccos(np.dot(n, e) / (np.linalg.norm(n)*e_mag)) * 180 / np.pi
            if e[2] < 0: arg_periapse = 360 - arg_periapse    # quadrant check


        # Calculate longitude of ascending node 
        # (angle between I and line of nodes)
        if i == 0: long_AN = None
        else:
            long_AN = np.arccos(n[0]/np.linalg.norm(n)) * 180 / (np.pi)
            if n[1] < 0: long_AN = 360 - long_AN   # quadrant check

        # Calculate the orbit period
        if energy >= 0:
            tp = None
        else:
            tp = 2*np.pi/np.sqrt(self.cb['mu']) * a**(3/2)

        # Calculate true anamoly of initial condiitons 
        # (angle between periapse and current position)
        if e_mag == 0: true_anamoly = None
        else:
            true_anamoly = np.arccos(np.dot(e, self.r0) / (np.linalg.norm(e)*np.linalg.norm(self.r0))) * 180 / np.pi
            if np.dot(self.r0, self.v0) < 0: true_anamoly = 360 - true_anamoly  # quadrant check

        # Calculate true longitude of periapse
        # (angle between I and eccentricity)
        if e_mag == 0: true_long_periapse = None
        else:
            true_long_periapse = np.arccos(e[0]/np.linalg.norm(e)) * 180 / np.pi
            if e[1] < 0: true_long_periapse = 360 - true_long_periapse     # quadrant check

        # Calculate the argument of latitude at epoch
        # (angle between line of nodes and satellite position)
        if i == 0: u = None
        else:
            u = np.arccos(np.dot(n, self.r0) / (np.linalg.norm(n)*np.linalg.norm(self.r0))) * 180 / np.pi
            if self.r0[2] < 0: u = 360 - u  # quadrant check

        # Calculate the true longitude at epoch
        # (angle between I and satellite position)
        true_long_epoch = np.arccos(self.r0[0] / np.linalg.norm(self.r0)) * 180 / np.pi
        if self.r0[1] < 0: true_long_epoch = 360 - true_long_epoch  # quadrant check

        if print_val:
            print('a = ' + str(a))
            print('e = ' + str(e_mag))
            print('i = ' + str(i))
            print('\u03C9 = ' + str(arg_periapse))
            print('\u03A9 = ' + str(long_AN))
            print('period = ' + str(tp))
            print('\u03BD = ' + str(true_anamoly))
            print('\u03A0 = ' + str(true_long_periapse))
            print('u = ' + str(u))
            print('l = ' + str(true_long_epoch))


        return np.array([a, e_mag, i, arg_periapse, long_AN, tp, true_anamoly, true_long_periapse, u, true_long_epoch])