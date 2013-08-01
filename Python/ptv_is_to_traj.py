#!/usr/bin/python
""" 
PTV_IS_TO_TRAJ converts a list of ptv_is.* files in the DIRECTORY into 
structure TRAJECTORY (shortcut is TRAJ). It's a simple bookkeeping procedure,
translated from Matlab ptv_is_to_traj.m
"""

import os
import glob
import numpy as np
# from numpy.lib.recfunctions import drop_fields,append_fields
from scipy.interpolate import splprep, splev
# import math as m
# import pdb


# --------------
def _read_ptv_is_files(directory,first=None,last=None):
	""" 
	_read_ptv_is_files(directory='.',first=None,last=None) 
	reads ptv_is.* files from the directory, ptv_is.first to ptv_is.last
	Example:
	data = _read_ptv_is_files('~/Desktop/GUI/pyptv2/test_for_v1_02/res',10000,100010)
	"""
	d = glob.glob(os.path.join(directory,'ptv_is.*'))
	d.sort()
	dlist = [int(os.path.splitext(name)[1][1:]) for name in d]
	if first is None or first < dlist[0]:
		first = 0
	else:
		first = dlist.index(first)
#
	if last is None or last > dlist[-1]:
		last = -1
	else:
		last = dlist.index(last)
	
	
	d = d[first:last] # only files between first and last
	dlist = dlist[first:last]
#	data = np.recarray(len(dlist),names='prev,next,x,y,z',formats="i4,i4,f8,f8,f8")
	data = []
	
	for ind,file in enumerate(d):
		print file
		frame = np.atleast_1d(np.genfromtxt(file,dtype="i4,i4,f8,f8,f8",
		names=['prev','next','x','y','z'],skip_header=1)).view(np.recarray)
		# print (ind, file, frame)
		frame.t = dlist[ind]*np.ones_like(frame.x)
		frame.trajid = np.atleast_1d(np.zeros_like(frame.x,dtype='i4'))
		data.append(frame)
	
	return data
	

#------------------------------
def _build_trajectories(data):
	"""
	build_trajectories(data) is responsible for the book keeping of the trajectories, 
	using the prev,next fields in the frames,creating new set of data, in the form of trajectories
	"""
	# first frame, initialization
	trajid = 0 # running trajid counter
	ind, = (data[0].next != -2).nonzero()
	if len(ind) > 0:
		data[0].trajid[ind] = range(trajid,trajid+len(ind))
		trajid += len(ind)
	for i,frame in enumerate(data[1:]):
		old, = (frame.prev > -1).nonzero() # those that continue from the previous frame
		frame.trajid[old] = data[i].trajid[frame.prev[old]-1] # notice i runs from zero
		# new trajectories - no background, but there's future
		ind, = np.logical_and(frame.prev == -1, frame.next != -2).nonzero()
		if len(ind) > 0:
			frame.trajid[ind] = range(trajid,trajid+len(ind))
			trajid += len(ind)
		frame.t = frame.t*np.ones_like(frame.x)
		
	return(data)
#----------------------------
def _ptv_to_traj(data,minlength=0):
	""" 
	ptv_to_traj(data) converts the data in frames with already prepared trajid to the data in form of
	trajectories
	"""
	# spline parameters
	s = 3.0 # smoothness parameter
	k = 3 # spline order
	nest = -1 # estimate of number of knots needed (-1 = maximal)
	
	
	
	trajid = np.hstack([f.trajid for f in data])
	x = np.hstack([f.x for f in data])
	y = np.hstack([f.y for f in data])
	z = np.hstack([f.z for f in data])
	t = np.hstack([f.t for f in data])

	traj = []
	ind, = ((x != 0) & ( y != 0) & (z != 0)).nonzero()

	x = x[ind]
	y = y[ind]
	z = z[ind]
	t = t[ind]
	trajid = trajid[ind]
	# import pdb; pdb.set_trace()
	dt = np.diff(np.unique(t)).min()
	unique_trajid = np.unique(trajid)
	print unique_trajid

	counter = -1
	trajlen = []
	for tr in unique_trajid:
		id = trajid == tr
		if sum(id) > minlength:
			counter += 1
			trajlen.append(sum(id))
			# find the knot points
			tckp,tmp = splprep([x[id],y[id],z[id]],s=s,k=k,nest=-1)
			new = np.linspace(0,1,trajlen[counter])
			# evaluate spline, including interpolated points
			xnew,ynew,znew = splev(new,tckp)
			u,v,w = splev(new,tckp,der=1)
			# 
			# ax,ay,az = splev(new,tckp,der=2)
			ax = np.gradient(u.flatten())
			ay = np.gradient(v.flatten())
			az = np.gradient(w.flatten())
# 			pdb.set_trace()
			tmp = np.rec.fromarrays([xnew,ynew,znew,u,v,w,ax,ay,az,t[id]],names='x,y,z,u,v,w,ax,ay,az,t')
			tmp.trajid = tr
			traj.append(tmp)

	return traj

# -------------------------
if __name__ == '__main__':
	print(' Reading ... \n')
	# directory = '/Volumes/alex/openptv_experiment/res';
	directory = '/Users/alex/Desktop/crowd_tracking/test_folder/res/'
	data = _read_ptv_is_files(directory,45160,45179) #('/Users/alex/Desktop/GUI/pyptv2/test_for_v1_02/res/')
	np.save('data',data)
	print(' Building trajectories ... \n')
	data = _build_trajectories(data)
	np.save('data',data)
	print('Restructuring and deriving velocity ... \n')
	traj = _ptv_to_traj(data)
	print(' Saving ... \n')
	np.save('traj',traj)
# 	
# 	% If manual cleaning is needed
# 	traj = graphical_cleaning_traj(traj,'xy')
# 	
# 	plot_long_trajectories(traj,10)
# 	
# 	% save(fullfile(directory,'traj.mat'),'traj')
# 	% load(fullfile(directory,'traj.mat'),'traj')
# 	
# 	s = 100;
# 	for i = 1:length(traj)
# 		newtraj(i) = link_trajectories_smoothn(traj(i),s);
# 	end
# 	
# 	% for i = 1:length(traj)
# 	%     newtraj(i) = link_trajectories_rbf(traj(i),0.1);
# 	% end
# 	[traj,removelist] = haitao_linking_criteria(newtraj,10,1.5);
# 	% save colloid_traj traj
# 	
# 	
# 	plot_long_trajectories(traj,5)
# 	hold on
# 	plot_long_trajectories(newtraj(removelist),1,1)