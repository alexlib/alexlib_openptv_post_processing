"""
Read ORI files
"""
import numpy as np

def read_ori(filename):
	f = open(filename,'r')
	# two first rows are camera positions: x,y,z and angles
	x,y,z = np.double(f.readline().split())
	print x,y,z
	alpha,beta,gamma = np.double(f.readline().split())
	print alpha,beta,gamma
	
	f.readline() # skip empty line
	
	# read transformation matrix
	mat = np.empty((3,3))
	print mat
	for i in range(3):
		mat[i,:] = np.array(np.double(f.readline().split()))
	
	print mat
	
	f.readline() # skip empty row
	
	# xp,yp are imaging center of the camera
	xp,yp = np.double(f.readline().split())
	print xp,yp
	
	# focal is the back-focal distance 
	focal = np.double(f.readline().split())
	print focal
	
	f.readline() # skip
	
	# vector pointing to the interface (glass window)
	
	try:
		xw,yw,zw = np.double(f.readline().split())
		print xw,yw,zw
	except:
		xw=yw=zw=[]
		print "old version"


	
	# still not clear what to return
	
	return x,y,z,alpha,beta,gamma,mat,xp,yp,focal,xw,yw,zw
	
