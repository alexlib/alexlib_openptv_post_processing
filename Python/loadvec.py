#!/usr/bin/pythonw
""" LOADVEC 

"""

import numpy as np
import matplotlib.pylab as pl


def get_header(fname):
	fname = 'Re1960/m002500.T000.D000.P001.H002.L.vec'
	f = open(fname)
	header = f.readline()
	f.close()
	ind1 = header.find('MicrosecondsPerDeltaT')
	dt = float(header[ind1:].split('"')[1])
	return dt

def get_data(fname):
	data = np.genfromtxt(fname,skip_header=1,delimiter=',',usecols=(0,1,2,3))
	
def read_directory(dirname):
	list_files = os.listdir(dirname)
	return list_files
	
	
