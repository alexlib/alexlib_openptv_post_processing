#!/usr/bin/pythonw
""" 
plot_long_trajectory(matfile) 
	: loads and plots the long trajectories from the LongTrajectory MAT file matfile
"""
from scipy.io import loadmat
import numpy as np
import pylab as pl
from mayavi import mlab
import sys
import os


filename = sys.argv[1]

data = loadmat(filename)
traj = data['traj'] # structure with fields: xf,yf,...

# @mlab.show
for d in traj:
	if d['xf'][0].size > 10:
		xf = d['xf'][0]
		yf = d['yf'][0]
		zf = d['zf'][0]
		u = (d['uf'][0]**2 + d['vf'][0]**2 + d['wf'][0]**2)**0.5		
		u[np.isinf(u)] = np.nan
		if np.all(np.isnan(u)):
			uf = np.gradient(xf.flatten())
			vf = np.gradient(yf.flatten())
			wf = np.gradient(zf.flatten())
			u = (uf**2 + vf**2 + wf**2)**0.5
			u[np.isinf(u)] = np.nan

		mlab.points3d(xf,yf,zf,np.atleast_2d(u).T,colormap="jet",scale_factor=10)

mlab.colorbar()
# arr = mlab.screenshot()
# 
# pl.imshow(arr)
# pl.axis('off')
# pl.show()
