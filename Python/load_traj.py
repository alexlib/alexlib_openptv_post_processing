#!/usr/bin/python
"""
Load traj structure from Matlab MAT file
"""
import numpy as np
from scipy.io import loadmat
import matplotlib as mp
import pylab as plt
from mayavi import mlab

path = '/Users/alex/Desktop/dikla_scene63/'
filename = 'scene63_traj.mat'
filename = ''.join([path,filename])

d = loadmat(filename,struct_as_record=True)
traj = d['traj'][0]
num_traj ,  = traj.shape
print("num of trajectories =  %d" % num_traj)

traj_len = np.empty((num_traj,1))
for i in range(num_traj):
	traj_len[i] = traj[i]['xf'].size
	
	
# plt.hist(traj_len), plt.show()

"""
xf = []
yf = []
zf = []
uf = []

for i in range(num_traj):
	if traj_len[i] > 5:
		xf.append(traj[i]['xf'].tolist())
		yf.append(traj[i]['yf'].tolist())
		zf.append(traj[i]['zf'].tolist())
		
"""
xf,yf,zf,axf = []

for i in range(num_traj):
	if traj_len[i] > 15:
		xf = np.vstack([xf,traj[i]['xf']])
		yf = np.vstack([yf,traj[i]['yf']])
		zf = np.vstack([zf,traj[i]['zf']])
		axf = np.vstack([uf,traj[i]['axf']])
		# mlab.plot3d(traj[i]['xf'],traj[i]['yf'],traj[i]['zf'],traj[i]['axf'])
		
		

mlab.points3d(xf,yf,zf,scale_factor=0.7)
		
mlab.show()
		




