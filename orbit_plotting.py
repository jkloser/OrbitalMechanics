import numpy as np
import central_bodies
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

def plot(r):

earth_rad = 6378.0;
earth_mu = 398600.0;

def diff_eq(t,y,mu):
	rx, ry, rz, vx, vy, vz
	r = np.array([rx,ry,rz])

	norm_r = np.linalg.norm(r)

	ax, ay, az = -r*mu/norm_r**3

	return [vx, vy, vz, rx, ry, rz]


if __name__ == '__main__':
	r_mag = earth_rad + 500
	v_mag = npsqrt(earth_mu/r_mag)

	r0=[r_mag, 0, 0]
	v0=[0,v_mag,0]

	tspan = 100*60.0
	dt = 100.0
	n_steps = int(np.ceil(tspan/dt))

	ys = np.zeros(n_steps,6)
	ts = np.zeros(n_steps,6)

	y0 = r0 + v0
	y[0] = np.array(y0)
	step = 1

	solver=ode(diff_eq)
	solver.set_integrator('lsoda')
	solver.set_initial_value(y0, 0)
	solver.set_f_params(earth_mu)


	# propagate
	while solver.successful() and step<n_steps:
		solver.integrate(solver.t+dt)
		ts[step]=solver.t
		ys[step]=solver.y
		step+=1

	rs=ys[:,:3]
	plot(rs)