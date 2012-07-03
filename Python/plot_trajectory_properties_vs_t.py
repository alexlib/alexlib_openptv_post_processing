#!/usr/bin/pythonw
""" 
plot_trajectory_properties_vs_t(matfile) 
	: loads and plots the long trajectories from the LongTrajectory MAT file matfile
"""
from scipy.io import loadmat
import numpy as np
import pylab as pl
# from mayavi import mlab
import sys
import os
from xuap_to_traj import Traj
from load_traj_mat import load_traj_mat

	
minLength = 6
eps_disp = 0.2


def convert_mat_to_traj(data):
	g = lambda name: name+" = row[0]['"+name+"'][0][0]"
	traj = []
	names = data[0][0].dtype.names 
	for row in data:
		for name in names:
			exec(g(name))
	
		traj.append(Traj(x,y,z,u,v,w,ax,ay,az,t,trajid))
	
	return traj



def plot_trajectory_vs_t(traj):
	t0 = traj[0].t[0]
	print t0
	print len(traj)
	
	pl.figure()
	ax1 = pl.subplot(311)
	pl.hold(True)
	ax2 = pl.subplot(312)
	pl.hold(True)
	ax3 = pl.subplot(313)
	pl.hold(True)
	
	for tr in traj:
		if tr.x.size > minLength:
# 			if (tr.x[-1] - tr.x[0])**2 + (tr.y[-1]-tr.y[0])**2 + (tr.z[-1]-tr.z[0])**2 < eps_disp:
# 				break
			
			t = tr.t - t0
			ax1.plot(tr.t,tr.x,'r-',tr.t,tr.y,'g--',tr.t,tr.z,'b-.')
			ax2.plot(tr.t,tr.u,'r-',tr.t,tr.v,'g--',tr.t,tr.w,'b-.')
			ax3.plot(tr.t,tr.ax,'r-',tr.t,tr.ay,'g--',tr.t,tr.az,'b-.')
			
	
	pl.show()


if __name__ == "__main__":
	if len(sys.argv) < 2:
		filename = 'marksingleparticle.mat'
	else:
		filename = sys.argv[1]

	
	file_created_with_Matlab = 0
	
	if file_created_with_Matlab:
		data = loadmat(filename)
		traj = data['traj'] # structure with fields: xf,yf,...
	else:
		traj = load_traj_mat(filename)
	

	plot_trajectory_vs_t(traj)
