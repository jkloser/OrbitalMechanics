import numpy as np
### Double check which units angles use, need to convert calculations to degrees

class orbit:
    def __init__(self, mu):
        self.mu = mu
        self.r0 = np.empty(3)
        self.v0 = np.empty(3)
        self.energy = np.empty(1)
        self.h = np.empty(3)
        self.e = np.empty(3)

    def init_cond(self, r, v):
        self.r0 = r
        self.v0 = v
        self.energy = np.linalg.norm(self.v0)**2/2 - self.mu/np.linalg.norm(self.r0)
        self.h = np.cross(self.r0, self.v0)
        self.e = np.cross(self.v0, self.h)/self.mu - self.r0/np.linalg.norm(self.r0)
        self.n = np.array([-h[2], h[1], 0])

    def calc_a(self):
        a = -self.mu/(2*self.energy)
        return a

    def calc_e(self):
        return np.linalg.norm(self.e)

    def calc_i(self):
        i = np.arccos(h[3]/np.linalg.norm(h))
        return i 

    def calc_arg_periapse(self):
        arg_periapse = np.arccos(np.dot(self.n, self.e) / (np.linalg.norm(self.n)*self.calc_e()))
        arg_periapse = arg_periapse * 360 / (2*np.pi)
        if self.e[3] < 0:
            arg_periapse = 360 - arg_periapse
        return arg_perigee

    def calc_long_AN(self):
        long_AN = np.arccos(self.n[1]/np.linalg.norm(self.n))
        long_AN = long_AN * 360 / (2*np.pi)
        if self.n[2] < 0:
            long_AN = 360 - long_AN
        return long_AN

    def calc_time_period(self):
        tp = 2*np.pi/np.sqrt(self.mu) * (self.mu/(2*self.energy))**(3/2)
        return tp

    ### Add funtionality such that if no new radius or vector is provided, the initial radius and vector are used.

    def calc_true_anamoly0(self):
        true_anamoly = np.arccos(np.dot(self.e, self.r0) / (np.linalg.norm(self.e)*np.linalg.norm(self.r0)))
        true_anamoly = true_anamoly * 360 / (2*np.pi)
        if np.dot(self.r0, self.v0) < 0:
            true_anamoly = 360 - true_anamoly
        return true_anamoly

    def calc_true_long_periapse(self):
        true_long_periapse = self.e[1]/np.linalg.norm(self.e)
        true_long_periapse = true_long_periapse * 360 / (2*np.pi)
        if self.e[2] < 0:
            true_long_periapse = 360 - true_long_periapse
        return true_long_periapse

    def calc_arg_lat_epoch(self):
        u = np.arccos(np.dot(self.n, self.r0) / (np.linalg.norm(self.n)*np.linalg.norm(self.r0)))
        u = u * 360 / (2*np.pi)
        if self.r0[3] < 0:
            u = 360 - u
        return u

    def calc_true_long_epoch(self):
        lambda = np.arccos(self.r0[1] / self.r0)
        lambda = lambda * 360 / (2*np.pi)
        if self.r0[2] < 0:
            lambda = 360 - lambda
        return lambda 