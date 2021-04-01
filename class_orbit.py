"""
John Kloser
AENG 4150 Orbital Mechanics
class_orbit.py
Project 1
"""

import numpy as np
import central_bodies
import matplotlib.pyplot as plt
#from scipy.integrate import ode
#from mpl_toolkits.mplot3d import Axes3D

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
        if round(i,self.tol) == 0 or round(e_mag,self.tol) == 0 or round(i-np.pi,self.tol) == 0: arg_periapse = float('NaN')
        else:
            arg_periapse = np.arccos(np.dot(n, e) / (np.linalg.norm(n)*e_mag))
            if e[2] < 0: arg_periapse = (2*np.pi) - arg_periapse    # quadrant check


        # Calculate longitude of ascending node 
        # (angle between I and line of nodes)
        if round(i,self.tol) == 0 or round(i-np.pi,self.tol) == 0: long_AN = float('NaN')
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
        if round(i,self.tol) == 0 or round(i-np.pi,self.tol) == 0: u = float('NaN')
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
            E = 2*np.arctanh(np.sqrt((e_mag-1)/(e_mag+1)) * np.tan(true_anamoly/2))
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

    def univ_formulation(self, dt, print_val=True, show_plot=True):
        r0_mag = np.linalg.norm(self.r0)
        v0_mag = np.linalg.norm(self.v0)
        energy = v0_mag**2/2 - self.cb['mu']/r0_mag
        a = -self.cb['mu']/(2*energy)

        alpha = (2*self.cb['mu']/r0_mag - v0_mag**2)/self.cb['mu'] # = 1/a

        # Add for loop here for an array of dt, initialize r and v
        if alpha > 0: # ellipse
            chi = np.sqrt(self.cb['mu']) * dt * alpha
        else: # hyperbola
            chi = np.sign(dt) * np.sqrt(-a) * np.log(-2*self.cb['mu']*alpha*dt / (np.dot(self.r0, self.v0)+np.sign(dt)*np.sqrt(self.cb['mu']*a)*(1-r0_mag/a)))

        # Iterating until the change in universal variable is less than tolerance
        iter = 1
        while True:
            z = chi**2*alpha

            if z > 0: # ellipse
                c = (1-np.cos(np.sqrt(z)))/z
                s = (np.sqrt(z)-np.sin(np.sqrt(z)))/np.sqrt(z**3)
            else: # hyperbola
                c = (1-np.cosh(np.sqrt(-z)))/z
                s = (np.sinh(np.sqrt(-z))-np.sqrt(-z))/np.sqrt((-z)**3)

            r = chi**2*c + np.dot(self.r0, self.v0)*chi*(1-z*s)/np.sqrt(self.cb['mu']) + r0_mag*(1-z*c)
            dchi = (np.sqrt(self.cb['mu'])*dt - chi**3*s - np.dot(self.r0, self.v0)*chi**2*c/np.sqrt(self.cb['mu']) - r0_mag*chi*(1-z*s))/r

            if round(dchi, self.tol) == 0.0:
                chi = chi + dchi
                print('iterations: ' + str(iter))
                break
            else:
                chi = chi + dchi
                iter+=1

        # f and g functions
        f = 1-chi**2*c/r0_mag
        g = dt - chi**3*s/np.sqrt(self.cb['mu'])
        gdot = 1 - chi**2*c/r
        fdot = np.sqrt(self.cb['mu'])*chi*(z*s-1)/(r*r0_mag)

        # Check accuracy and compute new radius/velocity
        check = f*gdot-fdot*g - 1
        r1 = f*self.r0 + g*self.v0
        v1 = fdot*self.r0 + gdot*self.v0

        if print_val:
            if round(check, self.tol) != 0:
                print('f and g check not valid')

            print('final X = ' + str(chi))
            print('f = ' + str(f))
            print('g = ' + str(g))
            print('fdot = ' + str(fdot))
            print('gdot = ' + str(gdot))
            print('check = ' + str(check))
            print()
            print('r1 = [' + str(r1[0]) + ', ' + str(r1[1]) + ', ' + str(r1[2]) + ']')
            print('v1 = [' + str(v1[0]) + ', ' + str(v1[1]) + ', ' + str(v1[2]) + ']')

        if show_plot:
            self.plot_PQW(r1, dt)

        return r1, v1

    def rot_PQW2IJK(self):
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

    def plot_PQW(self, points=np.array([]), dt=[], elements = []):
        # Get orbital elements if not provided
        if elements == []:
            elements = self.rv2elem(print_val=True)
        e_mag = elements[1]
        p = elements[0]*(1-e_mag**2)
        true_long_periapse = elements[7]

        # Create angles for calculating polar coords
        theta = np.linspace(0, 2*np.pi, 200)
        # Central body plot
        x_cb = self.cb['radius']*np.cos(theta)
        y_cb = self.cb['radius']*np.sin(theta)

        fig,ax = plt.subplots()
        ax.fill(x_cb, y_cb, "b", alpha=0.5)

        # Check if the true longitude of periapse (Pi) exists, if not, set to 0.
        # Indicates circular orbit
        if np.isnan(true_long_periapse):
            true_long_periapse = 0

        if e_mag > 1:
            # limit radius from going to infinity for hyperbolic orbit and reorder to prevent connecting two sides
            theta = theta[(e_mag*np.cos(theta)) >= -1]
            theta = np.concatenate((theta[theta>np.pi], theta[theta<np.pi])) # angles for continuous plot from bottom to top of hyperbola

        # Calculate the radii and the corresponding angles. Convert from polar to cartesian
        # True longitude at epoch accounts for location of periapsis wrt x-axis
        r_orbit = p/(1+e_mag*np.cos(theta))
        true_long_epoch = theta+true_long_periapse
        x_orbit = r_orbit*np.cos(true_long_epoch)
        y_orbit = r_orbit*np.sin(true_long_epoch)
        ax.plot(x_orbit, y_orbit, color='xkcd:lightish blue')

        # Add x- & y-axis arrows, eccentricity vector, annotations, axes units
        ax.arrow(0,0,0.5*abs(r_orbit[0]),0, head_width=0.1, head_length=0.1, fc='k', ec='k')
        ax.arrow(0,0,0,0.5*abs(r_orbit[0]), head_width=0.1, head_length=0.1, fc='k', ec='k')
        ax.annotate('x', xy=(0.55*abs(r_orbit[0]), 0), xycoords='data')
        ax.annotate('y', xy=(0, 0.55*abs(r_orbit[0])), xycoords='data')

        e_indx = np.where(r_orbit==min(r_orbit))
        ax.arrow(0,0,x_orbit[e_indx[0][0]],y_orbit[e_indx[0][0]], head_width = 0.02, head_length=0.02, fc='k', ec='k')
        ax.annotate('e', xy=(1.05*x_orbit[e_indx[0][0]],1.05*y_orbit[e_indx[0][0]]), xycoords='data')
        ax.set(xlabel=self.cb['units'], ylabel=self.cb['units'])

        if points.size > 1:
            points = np.vstack((self.r0, points))
            dt = np.vstack((0.00, dt))
            r_point = np.linalg.norm(points, axis=1)
            true_anamoly = np.arccos((p/r_point-1)/e_mag)

        else:
            r_point = np.linalg.norm(self.r0)
            true_anamoly = np.arccos((p/r_point-1)/e_mag)

        true_long_periapse_point = true_anamoly+true_long_periapse
        x_point = r_point*np.cos(true_long_periapse_point)
        y_point = r_point*np.sin(true_long_periapse_point)
        ax.plot(x_point, y_point, 'ro', alpha=0.7)

        for i, txt in enumerate(dt):
            ax.annotate(str(np.round(txt,2)), (x_point[i], y_point[i]))

        ax.axis('equal')
        plt.show()

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
    print(lst)
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

def gauss_problem(r1, r2, dt, DM, cb, tol = 5):
    r1_mag = np.linalg.norm(r1)
    r2_mag = np.linalg.norm(r2)
    delta_nu = np.dot(r1, r2) / (r1_mag*r2_mag)

    #A = np.sqrt(r1_mag*r2_mag)*np.sin(delta_nu) / np.sqrt(1-cos(delta_nu))

    DM = np.sign(np.pi-delta_nu)
    if DM>0: delta_nu = 2*np.pi - delta_nu
    A = DM*np.sqrt(r0_mag*r2_mag*(1+np.cos(delta_nu)))

    z = 0
    c = 1/2
    s = 1/6
    dt_iter = 0

    while round(dt-dt_iter, tol) != 0:
        y = r1_mag + r2_mag - A*(1-z*s)/np.sqrt(c)
        x = np.sqrt(y/c)
        dt_iter = (x**3*s + A*np.sqrt(y))/np.sqrt(cb['mu'])

        s_prime = 1/(2*z)*(c-3*s)
        c_prime = 1/(2*z)*(1-z*s-2*c)
        #if round(z, 3)==0:
        #    s_prime = 1/24
        #    c_prime = 1/120
        dtdz = (x**3 * (s_prime-(3*s*c_prime)/(2*c)) + A/8*(3*s*np.sqrt(y)/c + A/x)) / np.sqrt(cb['mu'])
        z = z + (dt-dt_iter)/dtdz

        c = 1/2 - z/np.math.factorial(4) + z**2/np.math.factorial(6) - x**3/np.math.factorial(8)
        s = 1/6 - z/np.math.factorial(5) + z**2/np.math.factorial(7) - x**3/np.math.factorial(9)

    f = 1-(y/r1_mag)
    g = A*np.sqrt(y/cb['mu'])
    gdot = 1-(y/r2_mag)

    v1 = (r2 - f*r1)/g
    v2 = (gdot*r2-r1)/g

    satellite = orbit(r1, v1, cb)
    satellite.rv2elem()

    return satellite, v2



