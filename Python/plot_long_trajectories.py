from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_long_trajectories(traj,minlength=10):
	""" Plots longer than minimum length trajectories in 3d
	Input:
	Output:
	Example:
	"""
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.hold(True)
	for d in traj:
		if len(d) > minlength:
			x = d.x
			y = d.y
			z = d.z
			ax.plot(x, y, z, label='parametric curve')
	
	plt.show()
	ax.hold(False)


if __name__ == '__main__':
	traj = sys.argv[1]
	plot_long_trajectories(traj)