import os
import glob
import numpy as np

def read_ptv_is_files(directory='.',first=None,last=None):
    fileList = glob.glob(directory+os.sep+'ptv_is.*')
    firstFile = int(fileList[0].split('ptv_is.')[-1])
    if first == None or first < firstFile:
        first = firstFile
	
    lastFile = int(fileList[-1].split('ptv_is.')[-1])
    if last == None or last > lastFile:
	last = lastFile
	
    dt = np.dtype([('prev',np.int),('next',np.int),('x',np.float),('y',np.float),('z',np.float)])
    tmp = np.loadtxt(fileList[0],skiprows=1,dtype=dt)
	
    numRows, = tmp.shape
    numFiles = lastFile - firstFile
	
    desc = np.dtype({'names' : ['data', 't'], 'formats': ['object', 'int']})
    traj = np.empty((numFiles,1),dtype=desc)
    traj[0] = np.array([(tmp,firstFile)],dtype=desc) 
	
    for k in range(1,numFiles):
        tmp = np.loadtxt(fileList[k],skiprows=1,dtype=dt)
        traj[k] = np.array([(tmp,firstFile+k)],dtype=desc) 
	
	return traj
	
	# traj[0]['data'][0]['next']
