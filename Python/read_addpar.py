"""
Read ADDPAR files
k1 [no unit], k2 [no unit], k3 [no unit], p1 [no unit], p2 [no unit], 
scale in x [no unit], shearing in [rad]
"""
import numpy as np

def read_addpar(filename):
	f = open(filename,'r')
	# two first rows are camera positions: x,y,z and angles
	try:
		k1,k2,k3,p1,p2,scale,shear = np.double(f.readline().split())
	except:
		k1=k2=k3=p1=p2=shear=0.0
		scale=1.0
		print "error reading addpar file, using defaults"

	return k1,k2,k3,p1,p2,scale,shear
	
