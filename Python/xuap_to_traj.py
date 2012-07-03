#!/usr/bin/python
import numpy as np
import glob
import os
from scipy.interpolate import splprep, splev
from mayavi import mlab
import scipy.io
import sys



class Frame:
	def __init__(self,prev,next,x,y,z,u,v,w,ax,ay,az,t,trajid=None):
		self.prev,self.next, self.x, self.y, self.z, self.u, self.v, self.w, self.ax, self.ay, \
		self.az, self.t, self.trajid = prev,next,x,y,z,u,v,w,ax,ay,az,t,trajid
	def length(self):
		return self.x.size
	
	def __call__(self):
		print(self.x,self.y,self.z)
#
	def save(self, filename=None):
		if filename is None:
			filename = 'temp.mat'
		scipy.io.savemat(filename, {'frames':self})
		
	def load(self, filename):
		scipy.io.loadmat(filename)


class Traj:
	def __init__(self,x,y,z,u,v,w,ax,ay,az,t,trajid):
		self.x, self.y, self.z, self.u, self.v, self.w, self.ax, self.ay, \
		self.az, self.t, self.trajid = x,y,z,u,v,w,ax,ay,az,t,trajid 
	
	def length(self):
		return self.x.size
		
	def plot3d(self):
		mlab.points3d(self.x,self.y,self.z,self.u**2+self.v**2+self.z**2,colormap="jet",scale_factor=10)
		
	def save(self, filename=None):
		if filename is None:
			filename = 'temp.mat'
		scipy.io.savemat(filename, {'traj':self})
		
	def load(self, filename):
		scipy.io.loadmat(filename)
	
	def __call__(self):
		return self.x, self.y, self.z, self.u, self.v, self.w, self.ax, self.ay, self.az, self.t, self.trajid



	

def read_xuap_files(directory,first=None,last=None):
# """ 
# read_xuap_files(directory='.',first=None,last=None) 
# reads xuap.* files from the directory, xuap.first to xuap.last
# Example:
#	data = _read_xuap_files('~/Desktop/GUI/pyptv2/test_for_v1_02/res',10000,100010)
# """
	d = glob.glob(os.path.join(directory,'xuap.*'))
	dlist = [int(os.path.splitext(name)[1][1:]) for name in d]
	if first is None or first < dlist[0]:
		first = dlist[0]
#
	if last is None or last > dlist[-1]:
		last = dlist[-1]
	
	
	d = d[dlist.index(first):dlist.index(last)] # only files between first and last
	dlist = dlist[dlist.index(first):dlist.index(last)]
	
#	print d
#	print dlist
#	data = n.recarray(len(dlist),names='prev,next,x,y,z',formats="i4,i4,f8,f8,f8")
	data = []
	
	for ind,f in enumerate(d):
#		print ind,f
		tmp = np.genfromtxt(f,dtype="i4,i4,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8",\
			names=['prev','next','xr','yr','zr','xf','yf','zf','u','v','w','ax','ay','az'],usecols=tuple(range(14))).view(np.recarray)
#		print frame
		frame = Frame(tmp.prev,tmp.next,tmp.xr,tmp.yr,tmp.zr,tmp.u,tmp.v,tmp.w,tmp.ax,tmp.ay,tmp.az,dlist[ind])
		data.append(frame)
	
	return data
# -----------------------------------
def build_trajectories_from_xuap(data):
	# first frame, initialization
	trajid = np.array([0],dtype='i4') # running trajid counter
	data[0].trajid = np.atleast_1d(np.nan*np.empty_like(data[0].x))
	data[0].t = np.atleast_1d(data[0].t*np.ones_like(data[0].x))
	ind = data[0].next > 0
	if np.atleast_1d(data[0].trajid).size > 1:
		data[0].trajid[ind.nonzero()] = range(trajid,trajid+ind.size)
	else:
		data[0].trajid = np.atleast_1d(trajid)
		
	trajid = trajid+ind.size
	for i,frame in enumerate(data[1:]):
		frame.trajid = np.atleast_1d(np.nan*np.empty_like(frame.x))
		old = np.atleast_1d(frame.prev) > 0
		if old.any():
			frame.trajid[old.nonzero()] = data[i].trajid[np.atleast_1d(frame.prev)[old.nonzero()]-1] # notice i runs from zero
		ind = np.logical_and(frame.prev == 0, frame.next > 0)
		if np.atleast_1d(frame.trajid).size > 1:
			frame.trajid[ind.nonzero()] = range(trajid,trajid+ind.size)
		else:
			frame.trajid = np.atleast_1d(trajid)
			
		trajid = trajid+ind.size
		frame.t = frame.t*np.ones_like(frame.x)
		
	return(data)
	
# -----------------
#----------------------------
def frame_to_traj(data,minlength=6):
	""" ptv_to_traj(data) converts the data in frames with already prepared trajid to the data in form of
	trajectories
	"""
	# spline parameters
	s=3.0 # smoothness parameter
	k=3.0 # spline order
	nest=-1 # estimate of number of knots needed (-1 = maximal)
	
	trajid = np.hstack([tmp.trajid.flatten() for tmp in data])
	x = np.hstack([tmp.x.flatten() for tmp in data])
	y = np.hstack([tmp.y.flatten() for tmp in data])
	z = np.hstack([tmp.z.flatten() for tmp in data])
	t = np.hstack([tmp.t.flatten() for tmp in data])
	ind = np.logical_or(np.logical_or(np.logical_or(x == 0.0, y == 0.0), z == 0.0),np.isnan(trajid))
	trajid = trajid[~ind].astype(np.int)
	x = x[~ind]
	y = y[~ind]
	z = z[~ind]
	t = t[~ind]

	dt = np.diff(np.unique(t)).min()
	traj = []
	counter = -1
	for tr in np.unique(trajid):
		id = trajid == tr
		trajlen = id.astype(np.int).sum()
		if  trajlen > minlength and np.diff(x[id].flatten()).any():
			counter += 1
			# find the knot points
			tckp,u = splprep([x[id],y[id],z[id]],s=s,k=k,nest=-1)
			new = np.linspace(0,1,trajlen)
			# evaluate spline, including interpolated points
			xnew,ynew,znew = splev(new,tckp)
			u,v,w = splev(new,tckp,der=1)
			# 
			# ax,ay,az = splev(new,tckp,der=2)
			ax = np.gradient(u.flatten())
			ay = np.gradient(v.flatten())
			az = np.gradient(w.flatten())
# 			print x,y,z,u,v,w,ax,ay,az,t[id],tr
			traj.append(Traj(xnew,ynew,znew,u,v,w,ax,ay,az,t[id],tr))

	return traj

if __name__ == "__main__":
	directory ='/Users/alex/Documents/Papers/Lid-Driven-Cavity/3dPTV_compression/Mark3DPTVSingleParticle/res'
	first =  104174 # 104235 #  #
	last =  146690 # 144335 # 
	data = read_xuap_files(directory,first,last)
	newdata = build_trajectories_from_xuap(data)
	traj = frame_to_traj(newdata)
	scipy.io.savemat('marksingleparticle.mat',{'traj':traj})
	
