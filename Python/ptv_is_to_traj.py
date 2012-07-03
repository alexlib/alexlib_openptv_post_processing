#!/usr/bin/python
""" 
PTV_IS_TO_TRAJ converts a list of ptv_is.* files in the DIRECTORY into 
structure TRAJECTORY (shortcut is TRAJ). It's a simple bookkeeping procedure,
translated from Matlab ptv_is_to_traj.m
"""

import os
import glob
import numpy as n
from numpy.lib.recfunctions import drop_fields,append_fields
from scipy.interpolate import splprep, splev
import math as m
# import pdb


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
		first = dlist[0]
	else:
		first = dlist.index(first)
#
	if last is None or last > dlist[-1]:
		last = dlist[-1]
	else:
		last = dlist.index(last)
	
	
	d = d[dlist.index(first):dlist.index(last)] # only files between first and last
	dlist = dlist[dlist.index(first):dlist.index(last)]
#	data = n.recarray(len(dlist),names='prev,next,x,y,z',formats="i4,i4,f8,f8,f8")
	data = []
	
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
	data[0].trajid = n.nan*n.empty_like(data[0].x)
	data[0].t = data[0].t*n.ones_like(data[0].x)
	ind = data[0].next > -2
	data[0].trajid[ind.nonzero()] = range(trajid,trajid+ind.size)
	trajid = trajid+ind.size
	for i,frame in enumerate(data[1:]):
		frame.trajid = n.nan*n.empty_like(frame.x)
		old = frame.prev > -1
		frame.trajid[old] = data[i].trajid[frame.prev[old]] # notice i runs from zero
		ind = n.logical_and(frame.prev < 0, frame.next > -2)
		frame.trajid[ind.nonzero()] = range(trajid,trajid+ind.size)
		trajid = trajid+ind.size
		frame.t = frame.t*n.ones_like(frame.x)
		
	return(data)
#----------------------------
def _ptv_to_traj(data,minlength=6):
# """ 
# ptv_to_traj(data) converts the data in frames with already prepared trajid to the data in form of
# trajectories
# """
	# spline parameters
	s=3.0 # smoothness parameter
	k=3.0 # spline order
	nest=-1 # estimate of number of knots needed (-1 = maximal)
	
	trajid = []
	trajlen = []
	x = []
	y = []
	z = []
	t = []
	for frame in data:
		for i,tr in enumerate(frame.trajid):
			if m.isnan(tr):
				pass
			else:
				trajid.append(n.int(tr))
				x.append(frame.x[i])
				y.append(frame.y[i])
				z.append(frame.z[i])
				t.append(frame.t[i])

	
	
	x = n.asarray(x)
	y = n.asarray(y)
	z = n.asarray(z)
	t = n.asarray(t)
	trajid = n.asarray(trajid)
# 	print x,y,z,trajid
	traj = []
	ind = n.logical_and(n.logical_and(x != 0.0, y != 0.0), z != 0.0)
	print ind
	x = x[ind]
	y = y[ind]
	z = z[ind]
	t = t[ind]
	trajid = trajid[ind]
	dt = n.diff(n.unique(t)).min()
	unique_trajid = n.unique(trajid)
	print unique_trajid
#traj = repmat(struct('xf',[],'yf',[],'zf',[],'uf',[],'vf',[],'wf',[],...
#    'axf',[],'ayf',[],'azf',[],'t',[],'trajid',[]),length(unique_trajid),1);
	counter = -1
	for tr in unique_trajid:
		id = trajid == tr
		if id.astype(n.int).sum() > minlength:
			counter += 1
			trajlen.append(id.astype(n.int).sum())
			# find the knot points
			tckp,u = splprep([x[id],y[id],z[id]],s=s,k=k,nest=-1)
			new = n.linspace(0,1,trajlen[counter])
			# evaluate spline, including interpolated points
			xnew,ynew,znew = splev(new,tckp)
			u,v,w = splev(new,tckp,der=1)
			# 
			# ax,ay,az = splev(new,tckp,der=2)
			ax = n.gradient(u.flatten())
			ay = n.gradient(v.flatten())
			az = n.gradient(w.flatten())
# 			pdb.set_trace()
			tmp = n.rec.fromarrays([xnew,ynew,znew,u,v,w,ax,ay,az,t[id]],names='x,y,z,u,v,w,ax,ay,az,t')
			tmp.trajid = tr
			traj.append(tmp)

	return traj

# -------------------------
if __name__ == '__main__':
	data = _read_ptv_is_files('/Users/alex/Desktop/res') #('/Users/alex/Desktop/GUI/pyptv2/test_for_v1_02/res/')
	data = _build_trajectories(data)
	traj = _ptv_to_traj(data)
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