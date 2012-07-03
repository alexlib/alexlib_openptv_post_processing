#!/usr/bin/python
import numpy as np
import glob
import os

def read_xuap_files(directory,first=None,last=None):
# """ 
# read_xuap_files(directory='.',first=None,last=None) 
# reads xuap.* files from the directory, xuap.first to xuap.last
# Example:
# 	data = _read_xuap_files('~/Desktop/GUI/pyptv2/test_for_v1_02/res',10000,100010)
# """
	d = glob.glob(os.path.join(directory,'xuap.*'))
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
	
# 	print d
# 	print dlist
#	data = n.recarray(len(dlist),names='prev,next,x,y,z',formats="i4,i4,f8,f8,f8")
	data = []
	
	for ind,f in enumerate(d):
# 		print ind,f
		frame = np.genfromtxt(f,dtype="i4,i4,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8",\
			names=['prev','next','xr','yr','zr','xf','yf','zf','u','v','w','ax','ay','az'],usecols=tuple(range(14))).view(np.recarray)
# 		print frame
		frame.t = dlist[ind]
		data.append(frame)
	
	return data
	
if __name__ == "__main__":
	directory = ='/Users/alex/Documents/Papers/Lid-Driven-Cavity/3dPTV_compression/Mark3DPTVSingleParticle/res'
	first = 104174
	last = 104184
	data = read_xuap_files(directory,first,last)
	np.savez('markxuapdata',data)