#!/usr/bin/python
""" load_traj_mat loads MAT files saved by xuap_to_traj.py and restores it as 
a list of Traj() class objects

Example:

>> traj = load_traj_mat('marksingleparticle.mat')


Author: Alex Liberzon
Date: 15-Dec-2010
"""
g = lambda name: name+" = row[0]['"+name+"'][0][0]"

from xuap_to_traj import Traj
from scipy.io import loadmat
import sys

 

def load_traj_mat(filename):
	data = loadmat(filename)['traj']
	traj = []
	names = data[0][0].dtype.names 
	for row in data:
		for name in names:
			exec(g(name))
	
		traj.append(Traj(x,y,z,u,v,w,ax,ay,az,t,trajid))
	
	return traj
	
if __name__ == "__main__":
	if len(sys.argv) < 2:
		filename = 'marksingleparticle.mat'
	else:
		filename = sys.argv[1]
		
	traj = convert_mat_to_traj(filename)