#!/usr/bin/python
""" 
PTV_IS_TO_TRAJ converts a list of ptv_is.* files in the DIRECTORY into 
structure TRAJECTORY (shortcut is TRAJ). It's a simple bookkeeping procedure,
translated from Matlab ptv_is_to_traj.m
"""

import os
import glob
import numpy as n
from numpy.lib.recfunctions import drop_fields
from scipy.interpolate import splprep, splev


# --------------
def _read_ptv_is_files(directory,first=None,last=None):
# """ 
# _read_ptv_is_files(directory='.',first=None,last=None) 
# reads ptv_is.* files from the directory, ptv_is.first to ptv_is.last
# Example:
# 	data = _read_ptv_is_files('~/Desktop/GUI/pyptv2/test_for_v1_02/res',10000,100010)
# """
	d = glob.glob(os.path.join(directory,'ptv_is.*'))
	dlist = [int(os.path.splitext(name)[1][1:]) for name in d]
	if first is None or first < dlist[0]:
		first = _first
	else:
		first = dlist.index(first)
#
	if last is None or last > dlist[-1]:
		last = _last
	else:
		last = dlist.index(last)
	
	data = []
	d = d[dlist.index(first):dlist.index(last)] # only files between first and last
	dlist = dlist[dlist.index(first):dlist.index(last)]
	
	for ind,file in enumerate(d):
		frame = n.genfromtxt(file,dtype="i4,i4,f8,f8,f8",names=['prev','next','x','y','z'],skip_header=1).view(n.recarray)
		frame.t = dlist[ind]
		data.append(frame)
	
	return data
	

#------------------------------
def _build_trajectories(data):
# """
# build_trajectories(data) is responsible for the book keeping of the trajectories, 
# using the prev,next fields in the frames,creating new set of data, in the form of trajectories
# """
	# first frame, initialization
	trajid = 0 # running trajid counter
	frame = data[0]
	frame.trajid = n.nan*n.empty_like(frame.x)
	ind = frame.next > -2
	frame.trajid[ind] = range(trajid,trajid+ind.size)
	trajid = trajid+ind.size


	for i,frame in data[1:]:
		frame.trajid = n.nan*n.empty_like(frame.x)
		old = frame.prev > -1
		frame.trajid[old] = data[i-1].trajid[frame.prev[old]]
		ind = frame.prev < 0 and frame.next > -2
		frame.trajid[ind] = range(trajid,trajid+ind.size)
		trajid = trajid+ind.size
		drop_fields(frame,['prev','next'])
		
	for frame in data:
		frame = frame[~n.isnan(frame)]
		frame.t = frame.t*n.ones_like(frame.x)
	
	return data
#----------------------------
def _ptv_to_traj(data,minlength=2):
# """ 
# ptv_to_traj(data) converts the data in frames with already prepared trajid to the data in form of
# trajectories
# """
	# spline parameters
	s=3.0 # smoothness parameter
	k=2 # spline order
	nest=-1 # estimate of number of knots needed (-1 = maximal)
	
	trajid = []
	x = []
	y = []
	z = []
	t = []
	for frame in data:
		trajid.extend(frame.trajid)
		x.extend(frame.x)
		y.extend(frame.y)
		z.extend(frame.z)
		t.extend(frame.t)
	
	x = n.asarray(x)
	y = n.asarray(y)
	z = n.asarray(z)
	t = n.asarray(t)
	trajid = n.asarray(trajid)
	traj = []
	ind = x != 0.0 and y != 0.0 and z != 0.0
	x = x[ind]
	y = y[ind]
	z = z[ind]
	t = t[ind]
	trajid = trajid[ind]
	dt = n.diff(n.unique(t)).min()
	unique_trajid = n.unique(trajid)
#traj = repmat(struct('xf',[],'yf',[],'zf',[],'uf',[],'vf',[],'wf',[],...
#    'axf',[],'ayf',[],'azf',[],'t',[],'trajid',[]),length(unique_trajid),1);
	counter = -1
	for k in unique_trajid:
		id = trajid == k
		if id.astype(n.int).sum() > minlength:
			counter += 1
			trajlen[counter] = id.astype(n.int).sum()
			# find the knot points
			tckp,u = splprep([x[id],y[id],z[id]],s=s,k=k,nest=-1)
			new = linspace(0,1,trajlen[counter])
			# evaluate spline, including interpolated points
			xnew,ynew,znew = splev(new,tckp)
			u,v,w = splev(new,tckp,der=1)
			ax,ay,az = splev(new,tckp,der=2)
			tmp = n.rec.fromarrays([xnew,ynew,znew,u,v,w,ax,ay,az,t[id]],names='x,y,z,u,v,w,ax,ay,az,t')
			tmp.trajid = k
			traj.append(tmp)

	return traj

# -------------------------
if __name__ == '__main__':

	data = _read_ptv_is_files('/Users/alex/Desktop/GUI/pyptv2/test_for_v1_02/res/')
	newdata = _build_trajectories(data)
	traj = _ptv_to_traj(newdata)
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