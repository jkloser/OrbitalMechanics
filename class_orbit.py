import numpy as np
import central_bodies
import matplotlib as plt
### Add a plot method
### Check fringe cases: undefined angles, parabola, hyperbola, circular orbits
### method to find r and v? given nu?

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
        # Calculate specific energy, specific angular momentum, eccentricity vector, line of nodes
        self.energy = np.linalg.norm(self.v0)**2/2 - self.cb['mu']/np.linalg.norm(self.r0)
        self.h = np.cross(self.r0, self.v0)
        self.e = np.cross(self.v0, self.h)/self.cb['mu'] - self.r0/np.linalg.norm(self.r0)
        self.n = np.array([-self.h[2], self.h[1], 0])

    def calc_a(self):
        # Calculate the semi-major axis. Returns 0 for a parabolic orbit.
        if self.energy == 0:
            a = 0
        else:
            a = -self.cb['mu']/(2*self.energy)
        return a

    def calc_e(self):
        # Calculate eccentricity magnitude
        return np.linalg.norm(self.e)

    def calc_i(self):
        # Calculate inclination of the orbit
        i = np.arccos(self.h[2]/np.linalg.norm(self.h)) * 180 / np.pi
        return i 

    def calc_arg_periapse(self):
        # Calculate argument of periapse
        # (angle between the line of nodes and eccentricity vector)
        arg_periapse = np.arccos(np.dot(self.n, self.e) / (np.linalg.norm(self.n)*self.calc_e()))
        arg_periapse = arg_periapse * 180 / (np.pi)
        
        # Quad check
        if self.e[2] < 0:
            arg_periapse = 360 - arg_periapse
        return arg_periapse

    def calc_long_AN(self):
        # Calculate longitude of ascending node 
        # (angle between I and line of nodes)
        long_AN = np.arccos(self.n[0]/np.linalg.norm(self.n))
        long_AN = long_AN * 180 / (np.pi)
        
        # Quad check
        if self.n[1] < 0:
            long_AN = 360 - long_AN
        return long_AN

    def calc_time_period(self):
        # Calculate the orbit period
        if self.energy == 0:
            return
        else:
            tp = 2*np.pi/np.sqrt(self.cb['mu']) * (self.cb['mu']/(2*self.energy))**(3/2)
            return tp

    def calc_true_anamoly(self):
        # Calculate true anamoly of initial condiitons 
        # (angle between periapse and current position)
        true_anamoly = np.arccos(np.dot(self.e, self.r0) / (np.linalg.norm(self.e)*np.linalg.norm(self.r0)))
        true_anamoly = true_anamoly * 180 / (np.pi)
        
        # Quad check
        if np.dot(self.r0, self.v0) < 0:
            true_anamoly = 360 - true_anamoly
        return true_anamoly

    def calc_true_long_periapse(self):
        # Calculate true longitude of periapse
        # (angle between I and eccentricity)
        true_long_periapse = self.e[0]/np.linalg.norm(self.e)
        true_long_periapse = true_long_periapse * 180 / (np.pi)

        # Quad check
        if self.e[1] < 0:
            true_long_periapse = 360 - true_long_periapse
        return true_long_periapse

    def calc_arg_lat_epoch(self):
        # Calculate the argument of latitude
        # (angle between line of nodes and satellite position)
        u = np.arccos(np.dot(self.n, self.r0) / (np.linalg.norm(self.n)*np.linalg.norm(self.r0)))
        u = u * 180 / (np.pi)

        # Quad check
        if self.r0[2] < 0:
            u = 360 - u
        return u

    def calc_true_long_epoch(self):
        # Calculate the true longitude at epoch
        # (angle between I and satellite position)
        true_long_epoch = np.arccos(self.r0[0] / np.linalg.norm(self.r0))
        true_long_epoch = true_long_epoch * 180 / (np.pi)

        # Quad check
        if self.r0[1] < 0:
            true_long_epoch = 360 - true_long_epoch
        return true_long_epoch