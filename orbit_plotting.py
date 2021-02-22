import numpy as np
import central_bodies
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

def plot(r):

earth_rad = 6378.0;
earth_mu = 398600.0;

def diff_eq(t,x,mu):

	"""
	 -   -        --
	|rxdot|      |rx|
	|rydot|      |ry|
	|rzdot|  = A |rz|
	|vxdot|      |vx|
	|vydot|      |vy|
	|vzdot|      |vz|
	 -   -        --
	 """

	# x is a state vector with the following elements
	r = np.array([x[0],x[1],x[2]])
	vx = x[3], vy = x[4], vz = x[5]

	norm_r = np.linalg.norm(r)
	ax, ay, az = -r*mu/norm_r**3

	# return state derivatives
	return [vx, vy, vz, ax, ay, az]


if __name__ == '__main__':
	r_mag = earth_rad + 500
	v_mag = npsqrt(earth_mu/r_mag)

	r0=[r_mag, 0, 0]
	v0=[0,v_mag,0]

	tspan = 100*60.0
	dt = 100.0
	n_steps = int(np.ceil(tspan/dt))

	states = np.zeros(n_steps,6)
	time = np.zeros(n_steps,6)

	init_state = np.concatenate(r0, v0)
	states[0] = init_state
	step = 1

	solver=ode(diff_eq)
	solver.set_integrator('lsoda')
	solver.set_initial_value(init_state, 0)
	solver.set_f_params(earth_mu)


	# propagate
	while solver.successful() and step<n_steps:
		solver.integrate(solver.t+dt)
		time[step]=solver.t
		state[step]=solver.x
		step+=1

	# rs=ys[:,:3]
	# plot(rs)