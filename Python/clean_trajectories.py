#!/usr/bin/pythonw
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

minLength = 10

filename = sys.argv[1]

data = loadmat(filename)
traj = data['traj'] # structure with fields: xf,yf,...
t0 = traj[0]['t'][0][0]

fig = pl.figure()
ax1 = pl.subplot(311)
pl.hold(True)
ax2 = pl.subplot(312)
pl.hold(True)
ax3 = pl.subplot(313)
pl.hold(True)

# @mlab.show
for d in traj:
	if d['xf'][0].size > minLength:
		xf = d['xf'][0]
		yf = d['yf'][0]
		zf = d['zf'][0]
		t =  d['t'][0] - t0
		uf = d['uf'][0]
		uf[np.isinf(uf)] = np.nan
		if np.all(np.isnan(uf)):
			uf = np.gradient(xf.flatten())
			vf = np.gradient(yf.flatten())
			wf = np.gradient(zf.flatten())
			axf = np.gradient(uf.flatten())
			ayf = np.gradient(vf.flatten())
			azf = np.gradient(wf.flatten())
		else:	
			vf = d['vf'][0]
			wf = d['wf'][0]
			axf = d['axf'][0]
			ayf = d['ayf'][0]
			azf = d['azf'][0]
				
		pl.axes(ax1)
		pl.plot(t,xf,'r-',t,yf,'g--',t,zf,'b-.')
		pl.axes(ax2)
		pl.plot(t,uf,'r-',t,vf,'g--',t,wf,'b-.')
		pl.axes(ax3)
		pl.plot(t,axf,'r-',t,ayf,'g--',t,azf,'b-.')
		

pl.show()
# arr = mlab.screenshot()
# 
# pl.imshow(arr)
# pl.axis('off')
# pl.show()
